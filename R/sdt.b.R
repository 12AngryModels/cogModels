
sdtClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "sdtClass",
    inherit = sdtBase,
    private = list(
        #### Init + run functions ----
        .run = function() {

            if (is.null(self$options$res) || is.null(self$options$stim))
                return()

            data <- private$.cleanData()
            results <- private$.compute(data)

            private$.preparePlot(results)
        },

        #### Compute results ----
        .compute = function(data) {

            path <- system.file("models/sdt.txt", package = "resWagner")

            # Sampling Parameters
            nChains = self$options$nChains
            nBurnin = self$options$nBurnin
            nSamples = self$options$nSamples
            nThin = self$options$nThin

            # Parameters to be monitored
            parameters <- c('dMu', 'cMu')

            # Draw posterior samples
            model <- rjags::jags.model(file=path, data=data, n.chains=nChains)

            if (nBurnin > 0)
                update(model, nBurnin)

            samples <- rjags::jags.samples(model, parameters, thin=nThin, n.iter=nSamples)

            r <- list(dMu=samples$dMu[,,], cMu=samples$cMu[,,])

            return(r)
        },

        #### Plot functions ----
        .preparePlot = function(results) {

            image <- self$results$plot

            density <- density(results$dMu)
            df <- data.frame(x = density$x, y = density$y)

            image$setState(df)

        },
        .plot = function(image, ggtheme, theme, ...) {

            if (is.null(image$state))
                return(FALSE)

            p <- ggplot2::ggplot(image$state, ggplot2::aes(x=x, y=y)) +
                ggplot2::geom_ribbon(ggplot2::aes(ymax=y), ymin=0, fill=theme$fill[2]) +
                ggplot2::geom_line(color=theme$color[1]) + ggplot2::xlab('dPrime') + ggplot2::ylab('Density') +
                ggtheme

            print(p)

            TRUE

        },

        #### Helper functions ----
        .cleanData = function() {

            df <- self$data

            response <- self$options$res
            stimulus <- self$options$stim
            signal <- self$options$sig

            subjects <- self$options$subj
            groups <- self$options$group

            nGroups <- 1
            if ( ! is.null(groups))
                nGroups <- length(unique(df[[groups]]))

            nSubjs <- 1
            if ( ! is.null(subjects)) {
                if (is.null(groups)) {
                    nSubjs <- length(unique(df[[subjects]]))
                } else {
                    nSubjs <- as.numeric(tapply(df[[subjects]],
                                                df[[groups]],
                                                function(x) return(length(unique(x)))))
                }
            }

            respData <- df[[response]]
            stimData <- df[[stimulus]]

            nSignal <- nNoise <- hits <- falseAlarms <- matrix(0, max(nSubjs), nGroups)

            dfs <- df
            if ( ! is.null(subjects)) {
                if ( ! is.null(groups)) {
                    dfs <- split(df, f = list(df[[subjects]], df[[groups]]))
                } else {
                    dfs <- split(df, f = list(df[[subjects]]))
                }
            }

            for (g in 1:nGroups) {
                for (i in 1:nSubjs[g]) {

                    sigTrials <- which(stimData == signal)
                    noiseTrials <- which(stimData != signal)

                    nSignal[i,g] <- length(sigTrials)
                    nNoise[i,g] <- length(noiseTrials)

                    hits[i,g] <- sum(respData[sigTrials] == stimData[sigTrials])
                    falseAlarms[i,g] <- sum(respData[noiseTrials] != stimData[noiseTrials])

                }
            }

            data <- list(nGroups = nGroups, nSubjs = nSubjs, nSignal = nSignal, nNoise = nNoise,
                         HR = hits, FAR = falseAlarms)

            print(data)

            return(data)
        })
)
