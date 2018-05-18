
sdtClass <- if (requireNamespace('jmvcore')) R6::R6Class(
    "sdtClass",
    inherit = sdtBase,
    private = list(
        .run = function() {
            
            path <- system.file("models/sdt.txt", package = "resWagner")
            
            res <- self$options$res
            trial <- self$options$trial
            
            if (is.null(res) || is.null(trial))
                return()
            
            private$.cleanData()
            
            # choiceData <- as.numeric(self$data[[choice]])
            # keyData <- as.numeric(self$data[[key]])
            # rewardData <- as.numeric(choiceData == keyData)
            # 
            # nTrials <- length(choiceData)
            # 
            # data <- list('choice'=choiceData, 'reward'=rewardData, 'N'=nTrials)
            # 
            # model <- rjags::jags.model(file = path,
            #                            data = data,
            #                            n.chains = 4,
            #                            n.adapt = 100)
            # 
            # update(model, 1000)
            # 
            # samples <- rjags::jags.samples(model, c('a'), 1000)
            # 
            # image <- self$results$plot
            # image$setState(samples[['a']])
            
        },
        .plot = function(image, ggtheme, theme, ...) {
            
            if (is.null(image$state))
                return(FALSE)
            
            plot(density(image$state))
            
            TRUE
            
        },
        .cleanData = function() {
            
            dataset <- self$data
            
            res <- self$options$res
            trial <- self$options$trial
            signal <- self$options$sig
            
            resData <- dataset[[res]]
            trialData <- dataset[[trial]]
            
            sigTrials <- which(trialData == signal)
            noiseTrials <- which(trialData != signal)
            
            HR <- sum(resData[sigTrials] == trialData[sigTrials])
            FAR <- sum(resData[noiseTrials] != trialData[noiseTrials])
            
            print(HR)
            
        })
)
