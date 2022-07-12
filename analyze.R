library(io)
library(ggplot2)

pheno <- qread("data/GDC-PANCAN.project_info.tsv");

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

samples.panels <- Reduce(
	intersect,
	lapply(counts.panel0s, function(counts.panel0) rownames(counts.panel0))
);

samples <- intersect(rownames(counts.exome0), samples.panels);

counts.exome <- counts.exome0[samples, ];
counts.panels <- lapply(counts.panel0s,
	function(counts.panel0) counts.panel0[samples, ]
);

opps.exome / sum(opps.exome)
lapply(opps.panels, function(opps.panel) opps.panel / sum(opps.panel))

tmc.exome <- rowSums(counts.exome);
tmb.exome <- tmc.exome / sum(opps.exome);
tmb.panels <- mapply(
	function(counts.panel, opps.panel) rowSums(counts.panel) / sum(opps.panel),
	counts.panels,
	opps.panels,
	SIMPLIFY=FALSE
);

par(mfrow=c(3, 2));
for (panel in panels) {
	plot(tmb.exome*1e6, tmb.panels[[panel]]*1e6, main=panel, log="xy");
	abline(a=0, b=1, col="red")
}

qdraw(
	{
		par(mfrow=c(3, 2));
		for (panel in panels) {
			smooth_scatter(tmb.exome*1e6, tmb.panels[[panel]]*1e6, main=panel);
			abline(a=0, b=1, col="red")
		}
	},
	width = 6,
	height = 9,
	file = "tmb-calib_tmb-exome-vs-panel.pdf"
);

qdraw(
	{
		par(mfrow=c(3, 2));
		for (panel in panels) {
			smooth_scatter(tmb.exome*1e6, tmb.panels[[panel]]*1e6, main=panel,
				xlim=c(0, 100), ylim=c(0, 100));
			abline(a=0, b=1, col="red")
		}
	},
	width = 6,
	height = 9,
	file = "tmb-calib_tmb-exome-vs-panel_zoom.pdf"
);

rmse <- function(x, y) {
	sqrt(mean((x - y)^2))
}

# @param x  reference
nrmse <- function(x, y) {
	rmse(x, y) / mean(x)	
}

fits.linear <- lapply(panels,
	function(panel) {
		lm(tmb.exome ~ tmb.panels[[panel]] - 1);
		#lm(tmb.exome ~ tmb.panels[[panel]]);
	}
);

tm.exposure <- rep(sum(opps.exome), length(tmc.exome));
fits.poisson <- lapply(panels,
	function(panel) {
		glm(tmc.exome ~ offset(log(tm.exposure)) + log(tmb.panels[[panel]]), family="quasipoisson");
	}
);

par(mfrow=c(3, 2));
for (panel in panels) {
	fit.linear <- fits.linear[[panel]];
	#plot(tmb.exome, fitted(fits.linear[[panel]]), main=panel)
	plot(tmb.exome, fitted(fits.linear[[panel]]), main=panel, xlim=c(0, 100/1e6), ylim=c(0, 100/1e6))
	points(tmb.exome, fitted(fits.poisson[[panel]]) / tm.exposure, col="blue")
	abline(a=0, b=1, col="red")
}

par(mfrow=c(3, 2));
for (panel in panels) {
	fit <- fits.linear[[panel]];
	smooth_scatter(tmb.exome, fitted(fit), main=panel, xlim=c(0, 100)/1e6, ylim=c(0, 100)/1e6)
	abline(a=0, b=1, col="red")
}

par(mfrow=c(3, 2));
for (panel in panels) {
	fit <- fits.poisson[[panel]];
	#plot(tmc.exome, fitted(fit), main=panel)
	smooth_scatter(tmc.exome, fitted(fit), main=panel, xlim=c(0, 5000), ylim=c(0, 5000))
	abline(a=0, b=1, col="red")
}

# ---

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

rates.exome <- lapply(mut.types,
	function(mut.type) {
		counts.exome[[mut.type]] / exposures.exome[[mut.type]]
	}
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

for (mut.type in mut.types) {
	par(mfrow=c(3, 2));
	for (panel in panels) {
		plot(rates.exome[[mut.type]], rates.panels[[panel]][[mut.type]], main=panel)
		abline(a=0, b=1, col="red")
	}
}

# mutation type specific Poisson model
fits.poisson.mut <- lapply(panels,
	function(panel) {
		lapply(mut.types,
			function(mut.type) {	
				exposures <- rep(exposures.exome[[mut.type]], length(counts.exome[[mut.type]]));
				glm(counts.exome[[mut.type]] ~
					offset(log(exposures)) + 
					log(rates.panels[[panel]][[mut.type]]),
					family="quasipoisson"
				)
			}
		)
	}
);

fits.poisson.mut.l2 <- lapply(panels,
	function(panel) {
		fits <- lapply(mut.types,
			function(mut.type) {	
				exposures <- rep(exposures.exome[[mut.type]], length(counts.exome[[mut.type]]));
				glm(counts.exome[[mut.type]] ~
					offset(log(exposures)) + 
					log(rates.panels[[panel]][[mut.type]]),
					family="quasipoisson"
				)
			}
		);

		# TODO replace with predictions on new data
		fitteds <- lapply(fits, fitted);

		# Include intercept term to allow second layer to fix
		# optimistic or pessimistic biases in rate estimates
		glm(
			#tmc.exome ~ c_a + c_g + c_t + t_a + t_c + t_g + indel - 1,
			tmc.exome ~ c_a + c_g + c_t + t_a + t_c + t_g + indel,
			data = fitteds,
			family="gaussian"
		);
	}
);

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

for (mut.type in mut.types) {
	par(mfrow=c(3, 2));
	for (panel in panels) {
		plot(counts.exome[[mut.type]], fitted(fits.poisson.mut[[panel]][[mut.type]]), main=panel)
		abline(a=0, b=1, col="red")
	}
}

fitted.pm <- function(fits) {
	preds <- lapply(mut.types, function(mut.type) fitted(fits[[mut.type]]));
	tmc.pred <- Reduce("+", preds);
}

par(mfrow=c(3, 2));
for (panel in panels) {
	#smooth_scatter(tmc.exome, fitted.pm(fits.poisson.mut[[panel]]), main=panel, xlim=c(0, 5000),  ylim=c(0, 5000))
	plot(tmc.exome, fitted.pm(fits.poisson.mut[[panel]]), main=panel)
	points(tmc.exome, fitted.pm(fits.poisson.mut.sep[[panel]]), main=panel, col="blue")
	points(tmc.exome, fitted(fits.poisson.mut.l2[[panel]]), main=panel, col="orange")
	abline(a=0, b=1, col="red")
}

tmc.exome.rep <- rep(tmc.exome, length(panels));
plot(tmc.exome.rep, fitted.pm(fits.poisson.mut.cpanel))
abline(a=0, b=1, col="red")

# ---

cd.linear <- unlist(lapply(fits.linear, function(fit) cor(tmb.exome, fitted(fit))^2));
cd.poisson <- unlist(lapply(fits.poisson, function(fit) cor(tmb.exome, fitted(fit)/tm.exposure)^2));
cd.poisson.mut <- unlist(lapply(fits.poisson.mut, function(fit) cor(tmb.exome, fitted.pm(fit)/tm.exposure)^2));
cd.poisson.mut.l2 <- unlist(lapply(fits.poisson.mut.l2, function(fit) cor(tmb.exome, fitted(fit)/tm.exposure)^2));
cd.poisson.mut.panelc <- cor(tmc.exome.rep, fitted.pm(fits.poisson.mut.cpanel))^2;

d <- rbind(
	data.frame(
		method = "linear",
		cd = cd.linear,
		panel = panels
	),
	data.frame(
		method = "poisson",
		cd = cd.poisson,
		panel = panels
	),
	data.frame(
		method = "poisson_mut",
		cd = cd.poisson.mut,
		panel = panels
	),
	data.frame(
		method = "poisson_mut_panelc",
		cd = cd.poisson.mut.panelc,
		panel = NA
	),
	data.frame(
		method = "poisson_mut_l2",
		cd = cd.poisson.mut.l2,
		panel = panels
	)
);

qdraw(
	ggplot(d, aes(x=method, y=cd)) + theme_classic() +
		geom_point(aes(shape=panel)) + coord_flip() +
		stat_summary(colour="red") +
		ylab("r^2") +
		theme(legend.position="bottom"),
	width = 6, height =4,
	file = "tmb-calib_training_r2.pdf"
)

nrmse.linear <- unlist(lapply(fits.linear, function(fit) nrmse(tmb.exome, fitted(fit))));
nrmse.poisson <- unlist(lapply(fits.poisson, function(fit) nrmse(tmb.exome, fitted(fit)/tm.exposure)));
nrmse.poisson.mut <- unlist(lapply(fits.poisson.mut, function(fit) nrmse(tmb.exome, fitted.pm(fit)/tm.exposure)));
nrmse.poisson.mut.l2 <- unlist(lapply(fits.poisson.mut.l2, function(fit) nrmse(tmb.exome, fitted(fit)/tm.exposure)));
nrmse.poisson.mut.panelc <- nrmse(tmc.exome.rep/tm.exposure, fitted.pm(fits.poisson.mut.cpanel)/tm.exposure);

d <- rbind(
	data.frame(
		method = "linear",
		nrmse = nrmse.linear,
		panel = panels
	),
	data.frame(
		method = "poisson",
		nrmse = nrmse.poisson,
		panel = panels
	),
	data.frame(
		method = "poisson_mut",
		nrmse = nrmse.poisson.mut,
		panel = panels
	),
	data.frame(
		method = "poisson_mut_panelc",
		nrmse = nrmse.poisson.mut.panelc,
		panel = NA
	),
	data.frame(
		method = "poisson_mut_l2",
		nrmse = nrmse.poisson.mut.l2,
		panel = panels
	)
);

qdraw(
	ggplot(d, aes(x=method, y=nrmse)) + theme_classic() +
		geom_point(aes(shape=panel)) + coord_flip() +
		stat_summary(colour="red") +
		theme(legend.position="bottom"),
	width = 6, height =4,
	file = "tmb-calib_training_nrmse.pdf"
)

