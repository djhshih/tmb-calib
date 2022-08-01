library(io)
library(ggplot2)

# --- Preamble

panels <- c(
	"MSK-IMPACT341", "MSK-IMPACT410", "MSK-IMPACT468", "MSK-IMPACT-HEME-399",
	"F1CDX", "DFCI"
);
names(panels) <- panels;

mut.types <- c("c_a", "c_g", "c_t", "t_a", "t_c", "t_g", "indel");
names(mut.types) <- mut.types;
counts.fields <- mut.types;
opps.fields <- c("a", "c", "g", "t");

smooth_scatter <- function(x, y, ...) {
	smoothScatter(x, y, transformation = function(x) x^0.15, ...)
}

rmse <- function(x, y) {
	sqrt(mean((x - y)^2))
}

# @param x  reference
nrmse <- function(x, y) {
	rmse(x, y) / mean(x)	
}

# --- Read input data

pheno0 <- qread("data/GDC-PANCAN.project_info.tsv");

counts.exome0 <- read.table("data/mut_count/Exome_nonsilent_mut.txt", sep="\t", row.names=1, header=FALSE);
colnames(counts.exome0) <- counts.fields;

counts.panel0s <- lapply(panels,
	function(panel) {
		counts.panel0 <- read.table(sprintf("data/mut_count/%s_nonsilent_mut.txt", panel), sep="\t", row.names=1, header=FALSE);
		colnames(counts.panel0) <- counts.fields;
		counts.panel0
	}
);

opps.exome <- read.table("data/nuc_content/Exome.final.bed.nuc", sep="\t", header=FALSE);
colnames(opps.exome) <- opps.fields;
opps.exome <- as.matrix(opps.exome)[1,];

opps.panels <- lapply(panels,
	function(panel) {
		opps.panel <- read.table(sprintf("data/nuc_content/%s.bed.nuc", panel), sep="\t", header=FALSE);
		colnames(opps.panel) <- opps.fields;
		as.matrix(opps.panel)[1,]
	}
);

# Ensure that samples are the same across platforms

samples.panels <- Reduce(
	intersect,
	lapply(counts.panel0s, function(counts.panel0) rownames(counts.panel0))
);

samples <- intersect(rownames(counts.exome0), samples.panels);

pheno <- pheno0[match(samples, pheno0$sample), ];
counts.exome <- counts.exome0[samples, ];
counts.panels <- lapply(counts.panel0s,
	function(counts.panel0) counts.panel0[samples, ]
);

# nucleotide distributions
opps.exome / sum(opps.exome)
lapply(opps.panels, function(opps.panel) opps.panel / sum(opps.panel))

# total mutation counts and burden
tmc.exome <- rowSums(counts.exome);
tmb.exome <- tmc.exome / sum(opps.exome);
tmb.panels <- mapply(
	function(counts.panel, opps.panel) rowSums(counts.panel) / sum(opps.panel),
	counts.panels,
	opps.panels,
	SIMPLIFY=FALSE
);

# Calculate the nucleotide-specific exposures

exposures.exome <- list(
	# C>A and G>T
	c_a = unname(opps.exome["c"] + opps.exome["g"]),
	# C>G and G>C
	c_g = unname(opps.exome["c"] + opps.exome["g"]),
	# C>T and G>A
	c_t = unname(opps.exome["c"] + opps.exome["g"]),
	# T>A and A>T
	t_a = unname(opps.exome["t"] + opps.exome["a"]),
	# T>C and A>G
	t_c = unname(opps.exome["t"] + opps.exome["a"]),
	# T>G and A>C
	t_g = unname(opps.exome["t"] + opps.exome["a"]),
	# indel
	indel = sum(opps.exome)
);

exposures.panels <- lapply(opps.panels,
	function(opps.panel) {
		cg <- opps.panel["c"] + opps.panel["g"];
		ta <- opps.panel["t"] + opps.panel["a"];
		list(
			c_a = cg, c_g = cg, c_t = cg,
			t_a = ta, t_c = ta, t_g = ta, indel = cg+ta
		)
	}
);

# Calculate type specific mutation rates

rates.exome <- lapply(mut.types,
	function(mut.type) {
		counts.exome[[mut.type]] / exposures.exome[[mut.type]]
	}
);

rates.panels <- mapply(
	function(counts.panel, exposures.panel) {
		mapply(
			function(x, n) {
				(x + 0.5) / (n + 1)
				# too little regularization (0.1, 0.2) => pessimistic bias
				# too much regularization (10, 20) => optimistic bias
			},
			counts.panel,
			exposures.panel,
			SIMPLIFY=FALSE
		)
	},
	counts.panels,
	exposures.panels,
	SIMPLIFY=FALSE
);

# Assess TMB estimated from each panel vs. exome

qdraw(
	{
		par(mfrow=c(3, 2));
		for (panel in panels) {
			smooth_scatter(tmb.exome*1e6, tmb.panels[[panel]]*1e6, main=panel,
				xlab="exome TMB per Mbp",
				ylab="panel TMB per Mbp"
			);
			abline(a=0, b=1, col="red")
		}
	},
	width = 6, height = 9,
	file = "tmb-calib_tmb-exome-vs-panel.pdf"
);

qdraw(
	{
		par(mfrow=c(3, 2));
		for (panel in panels) {
			smooth_scatter(tmb.exome*1e6, tmb.panels[[panel]]*1e6, main=panel,
				xlim=c(0, 100), ylim=c(0, 100),
				xlab="exome TMB per Mbp",
				ylab="panel TMB per Mbp"
			);
			abline(a=0, b=1, col="red")
		}
	},
	width = 6, height = 9,
	file = "tmb-calib_tmb-exome-vs-panel_zoom.pdf"
);

qdraw(
	{
		par(mfrow=c(3, 2));
		for (panel in panels) {
			smooth_scatter(log10(tmb.exome*1e6), log10(tmb.panels[[panel]]*1e6), main=panel,
				xlab="log10(exome TMB per Mbp)",
				ylab="log10(panel TMB per Mbp)"
			)
			abline(a=0, b=1, col="red")
		}
	},
	width = 6, height = 9,
	file = "tmb-calib_tmb-exome-vs-panel_log.pdf"
);

qdraw(
	{
		par(mfrow=c(6, 7));
		for (panel in panels) {
			for (mut.type in mut.types) {
				smooth_scatter(rates.exome[[mut.type]], rates.panels[[panel]][[mut.type]], main=paste(panel, mut.type))
				abline(a=0, b=1, col="red")
			}
		}
	},
	width = 16, height = 16,
	file = "tmb-calib_smb-exome-vs-panel.pdf"
);

# --- Models

train_linear <- function(tmb.exome, tmb.panel) {
	d <- data.frame(
		y = tmb.exome,
		x = tmb.panel
	);
	# no y-intercept because it can result in negative TMB predictions
	fit <- lm(y ~ x - 1, data=d);
	structure(list(model=fit), class="tmb_calib_linear")
}

# @param tmb.panel  TMB measured by a targeted panel
predict_linear <- function(model, tmb.panel) {
	d <- data.frame(
		x = tmb.panel
	);
	predict(model$model, d)
}

train_poisson <- function(tmc.exome, opps.exome, tmb.panel) {
	d <- data.frame(
		count = tmc.exome,
		log_exposure = log(sum(opps.exome)),
		log_x = log(tmb.panel)
	);
	fit <- glm(count ~ offset(log_exposure) + log_x, family="quasipoisson", data=d);

	params <- list(
		exposure = sum(opps.exome)
	);

	structure(list(model=fit, params=params), class="tmb_calib_poisson")
}

# @param tmb.panel  TMB measured by a targeted panel
predict_poisson <- function(model, tmb.panel) {
	d <- data.frame(
		log_exposure = log(model$params$exposure),
		log_x = log(tmb.panel)
	);
	# response is the count
	predict(model$model, d, type="response") / model$params$exposure
}

# mutation type specific Poisson model
train_poisson_smb <- function(counts.exome, exposures.exome, rates.panel) {
	mut.types <- names(rates.panel);
	names(mut.types) <- mut.types;

	fits <- lapply(
		mut.types,
		function(mut.type) {
			d <- data.frame(
				count = counts.exome[[mut.type]],
				log_exposure = log(exposures.exome[[mut.type]]),
				log_x = log(rates.panel[[mut.type]])
			);
			glm(count ~ offset(log_exposure) + log_x, family="quasipoisson", data=d)
		}
	);

	params <- list(
		exposure = exposures.exome
	);

	structure(list(models=fits, params=params), class="tmb_calib_poisson_smb")
}

# @param smb.panel  specific mutation burden measured by a targeted panel
predict_poisson_smb <- function(model, smb.panel) {
	preds <- lapply(
		names(smb.panel),
		function(mut.type) {
			d <- data.frame(
				log_exposure = log(model$params$exposure[[mut.type]]),
				log_x = log(smb.panel[[mut.type]])
			);
			# NB response is the count... why does predict not use offset?
			predict(model$models[[mut.type]], d, type="response")
		}
	);

	exposure <- model$params$exposure$indel;
	Reduce("+", preds) / exposure
}


# mutation type specific Poisson model
train_poisson_smb_2l <- function(counts.exome, exposures.exome, rates.panel) {
	fits <- train_poisson_smb(counts.exome, exposures.exome, rates.panel)$models;

	fitteds <- lapply(fits, fitted);
	tmc.exome <- rowSums(counts.exome);

	# include intercept term to allow second layer to fix
	# optimistic or pessimistic biases in rate estimates
	ensemble <- glm(
		tmc.exome ~ c_a + c_g + c_t + t_a + t_c + t_g + indel,
		data = fitteds,
		family="gaussian"
	);

	params <- list(
		exposure = exposures.exome
	);

	structure(list(models.l1=fits, model=ensemble, params=params),
		class="tmb_calib_poisson_smb_2l")
}

# @param smb.panel  specific mutation burden measured by a targeted panel
predict_poisson_smb_2l <- function(model, smb.panel) {
	mut.types <- names(smb.panel);
	names(mut.types) <- mut.types;

	preds <- lapply(
		mut.types,
		function(mut.type) {
			d <- data.frame(
				log_exposure = log(model$params$exposure[[mut.type]]),
				log_x = log(smb.panel[[mut.type]])
			);
			# NB response is the count
			predict(model$models.l1[[mut.type]], d, type="response")
		}
	);
	
	# NB response is the count
	# total exposure is the number of nucleotides in the reference panel (i.e.
	# exome), which is the same as the indel exposure
	exposure <- model$params$exposure$indel;
	pmax(0, predict(model$model, newdata=preds, type="response") / exposure)
}

# Obsolete models

# mutation type specific Poisson model 
# with separate terms for event counts and exposures
fits.poisson.mut.sep <- lapply(panels,
	function(panel) {
		lapply(mut.types,
			function(mut.type) {	
				n <- length(counts.exome[[mut.type]]);
				exposures <- rep(exposures.exome[[mut.type]], n);
				exposures.panel.mut <- rep(exposures.panels[[panel]][[mut.type]], n);
				glm(counts.exome[[mut.type]] ~
					offset(log(exposures)) + 
					log(counts.panels[[panel]][[mut.type]] + 0.5) +
					log(exposures.panel.mut + 1),
					family="quasipoisson"
				)
			}
		)
	}
);

fits.poisson.mut.cpanel <- lapply(mut.types,
	function(mut.type) {	
		N <- length(counts.exome[[mut.type]]);
		G <- length(panels);
		d <- data.frame(
			counts = rep(counts.exome[[mut.type]], G),
			exposure = rep(rep(exposures.exome[[mut.type]], N), G),
			panel = rep(unlist(panels), each=N),
			rate = unlist(lapply(panels, function(p) rates.panels[[p]][[mut.type]]))
		);
		glm(counts ~ offset(log(exposure)) + panel + log(rate) - 1,
			data=d, family="quasipoisson"
		)
	}
);

# --- Evaluate training performance

models.linear <- lapply(panels,
	function(panel) train_linear(tmb.exome, tmb.panels[[panel]])
);
tr.fitted.linear <- lapply(panels,
	function(panel) predict_linear(models.linear[[panel]], tmb.panels[[panel]])
);

models.poisson <- lapply(panels,
	function(panel) train_poisson(tmc.exome, opps.exome, tmb.panels[[panel]])
);
tr.fitted.poisson <- lapply(panels,
	function(panel) predict_poisson(models.poisson[[panel]], tmb.panels[[panel]])
);

models.poisson.smb <- lapply(panels,
	function(panel) train_poisson_smb(counts.exome, exposures.exome, rates.panels[[panel]])
);
tr.fitted.poisson.smb <- lapply(panels,
	function(panel) predict_poisson_smb(models.poisson.smb[[panel]], rates.panels[[panel]])
);

models.poisson.smb.2l <- lapply(panels,
	function(panel) train_poisson_smb_2l(counts.exome, exposures.exome, rates.panels[[panel]])
);
tr.fitted.poisson.smb.2l <- lapply(panels,
	function(panel) predict_poisson_smb_2l(models.poisson.smb.2l[[panel]],
		rates.panels[[panel]])
);

par(mfrow=c(3, 2));
for (panel in panels) {
	smooth_scatter(tmb.exome, tr.fitted.linear[[panel]], main=panel,
		xlim=c(0, 100)/1e6, ylim=c(0, 100)/1e6
	)
	abline(a=0, b=1, col="red")
}

par(mfrow=c(3, 2));
for (panel in panels) {
	smooth_scatter(tmb.exome, tr.fitted.poisson[[panel]], main=panel,
		xlim=c(0, 100)/1e6, ylim=c(0, 100)/1e6
	)
	abline(a=0, b=1, col="red")
}

par(mfrow=c(3, 2));
for (panel in panels) {
	smooth_scatter(tmb.exome, tr.fitted.poisson.smb[[panel]], main=panel,
		xlim=c(0, 100)/1e6, ylim=c(0, 100)/1e6
	)
	abline(a=0, b=1, col="red")
}

par(mfrow=c(3, 2));
for (panel in panels) {
	smooth_scatter(tmb.exome, tr.fitted.poisson.smb.2l[[panel]], main=panel,
		xlim=c(0, 100)/1e6, ylim=c(0, 100)/1e6
	)
	abline(a=0, b=1, col="red")
}

par(mfrow=c(3, 2));
for (panel in panels) {
	plot(tmb.exome, tr.fitted.linear[[panel]], main=panel,
		xlim=c(0, 100/1e6), ylim=c(0, 100/1e6)
	)
	points(tmb.exome, tr.fitted.poisson[[panel]], col="blue")
	points(tmb.exome, tr.fitted.poisson.smb.2l[[panel]], col="orange")
	abline(a=0, b=1, col="red")
}

qdraw(
	{
		par(mfrow=c(6, 7));
		for (panel in panels) {
			for (mut.type in mut.types) {
				smooth_scatter(
					counts.exome[[mut.type]],
					fitted(models.poisson.smb[[panel]]$models[[mut.type]]),
					main=paste(panel, mut.type)
				)
				abline(a=0, b=1, col="red")
			}
		}
	},
	width = 16, height = 16,
	file = "tmb-calib_poisson-smb_smc.pdf"
);

# @param yhatss  a list of list of fitted values
# @param f   evaluation function
evaluate_fits <- function(yhatss, f) {
	d <- do.call(rbind, mapply(
		function(yhats, method) {
			do.call(rbind, mapply(
				function(yhat, group) {
					data.frame(
						method = method,
						score = f(yhat),
						group = group
					)
				},
				yhats, names(yhats),
				SIMPLIFY = FALSE
			))
		},
		yhatss, names(yhatss),
		SIMPLIFY = FALSE
	));
	rownames(d) <- NULL;
	d$method <- factor(d$method, levels=unique(d$method));
	d$group <- factor(d$group, levels=unique(d$group));

	d
}

# ---

fitteds.tr <- list(
	none = tmb.panels,
	linear = tr.fitted.linear,
	poisson = tr.fitted.poisson,
	poisson_smb_2l = tr.fitted.poisson.smb.2l
);

d.tr.cd <- evaluate_fits(fitteds.tr, function(yhat) cor(tmb.exome, yhat)^2);

qdraw(
	ggplot(d.tr.cd, aes(x=method, y=score)) + theme_classic() +
		geom_point(aes(shape=group)) +
		geom_line(aes(group=group), alpha=0.2) +
		stat_summary(colour="red") +
		labs(shape="panel") +
		ylab("R^2") +
		theme(legend.position="bottom")
	,
	width = 6, height =4,
	file = "tmb-calib_training_r2.pdf"
)

d.tr.nrmse <- evaluate_fits(fitteds.tr, function(yhat) nrmse(tmb.exome, yhat));

qdraw(
	ggplot(d.tr.nrmse, aes(x=method, y=score)) + theme_classic() +
		geom_point(aes(shape=group)) +
		geom_line(aes(group=group), alpha=0.2) +
		stat_summary(colour="red") +
		ylab("NRMSE") +
		theme(legend.position="bottom")
	,
	width = 6, height =4,
	file = "tmb-calib_training_nrmse.pdf"
)

# --- Evaluate cross-validation performance

cv <- function(X, y, folds, train_f, predict_f) {
	if (is.vector(X)) {
		X <- matrix(X, ncol=1);
	}
	if (is.vector(y)) {
		y <- matrix(y, ncol=1);
	}
	fs <- unique(folds);
	yhat <- rep(NA, nrow(y));
	for (f in fs) {
			idx <- folds != f;
			fit <- train_f(X[idx, ], y[idx, ]);
			yhat[!idx] <- predict_f(fit, X[!idx, ])
	}
	yhat
}

#set.seed(1234);
#k <- 5;
#folds <- sample.int(k, size=nrow(pheno), replace=TRUE);

folds <- pheno$disease_code.project;

fitted.linear <- lapply(panels,
	function(panel) {
		cv(tmb.panels[[panel]], tmb.exome, folds,
			function(X, y) train_linear(y, X),
			function(fit, X.new) predict_linear(fit, X.new)
		);
	}
);

fitted.poisson <- lapply(panels,
	function(panel) {
		cv(tmb.panels[[panel]], tmc.exome, folds,
			function(X, y) train_poisson(y, opps.exome, X),
			function(fit, X.new) predict_poisson(fit, X.new)
		)
	}
);

fitted.poisson.smb.2l <- lapply(panels,
	function(panel) {
		cv(as.data.frame(rates.panels[[panel]]), counts.exome, folds,
			function(X, y) train_poisson_smb_2l(y, exposures.exome, X),
			function(fit, X.new) predict_poisson_smb_2l(fit, X.new)
		)
	}
);

# ---

fitteds.cv <- list(
	none = tmb.panels,
	linear = fitted.linear,
	poisson = fitted.poisson,
	poisson_smb_2l = fitted.poisson.smb.2l
);

d.cd <- evaluate_fits(fitteds.cv, function(yhat) cor(tmb.exome, yhat)^2);

qdraw(
	ggplot(d.cd, aes(x=method, y=score)) + theme_classic() +
		geom_point(aes(shape=group)) +
		geom_line(aes(group=group), alpha=0.2) +
		stat_summary(colour="red") +
		labs(shape="panel") +
		ylab("R^2") +
		theme(legend.position="bottom")
	,
	width = 6, height =4,
	file = "tmb-calib_cv_r2.pdf"
)

d.nrmse <- evaluate_fits(fitteds.cv, function(yhat) nrmse(tmb.exome, yhat));

qdraw(
	ggplot(d.nrmse, aes(x=method, y=score)) + theme_classic() +
		geom_point(aes(shape=group)) +
		geom_line(aes(group=group), alpha=0.2) +
		stat_summary(colour="red") +
		ylab("NRMSE") +
		theme(legend.position="bottom")
	,
	width = 6, height =4,
	file = "tmb-calib_cv_nrmse.pdf"
)

