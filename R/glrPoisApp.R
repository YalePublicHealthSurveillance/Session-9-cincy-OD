surv.ds.convert <- function(ds, datevar='date', casevar='cases'){
  year.start <- min(lubridate::year(ds[,datevar]))
  week.start <- min(lubridate::week(ds[,datevar])[lubridate::year(ds[,datevar])==year.start])
  
  SurvObj <- create.disProg(
    week = 1:nrow(ds), #Index of observations
    observed = ds[,casevar] ,  #cases
    state=matrix(0, nrow=nrow(ds), ncol=1), #just 0s
    start = c(year.start, week.start)) #start date; 1st week of 1990
  return(SurvObj)
}


glrpois_App <-function(ds=ds1, datevar='date', casevar='cases', n.weeks.train=53){
  surv.ds1 <- surv.ds.convert(ds, datevar=datevar, casevar=casevar)
  
  
  shinyApp(
    ui=fluidPage(
      
      
      sliderInput("week.test", "Current Week:",
                  min=(n.weeks.train), max=nrow(ds), value=(n.weeks.train), step=1),
      sliderInput("set.thresh", "Threshold for alarm (default=5):",
                  min=1, max=10, value=5),
      selectInput('adjust.season', 'Adjust seasonality?', selected='1 seasonal term', choices=c('No adjustment', '1 seasonal term', '2 seasonal terms')  ),
      checkboxInput('adjust.trend', 'Adjust trend?', value=F ),
      checkboxInput('fit.negbin', 'Use negative binomial instead of Poisson?', value=F ),
      
      plotOutput("periodPlot")
    ),
    server=function(input, output){
      output$periodPlot = renderPlot({
        
        season.adjust <- 0
        season.adjust[input$adjust.season=='1 seasonal term'] <- 1
        season.adjust[input$adjust.season=='2 seasonal terms'] <- 2
        
        if(input$fit.negbin==T){
          mod1<- algo.glrnb(surv.ds1,control=list(
            range=c(n.weeks.train:input$week.test),
            c.ARL=input$set.thresh,
            M=-1, #How many time points back should we look? Negative 1: use all cases
            ret=c('value'),
            alpha=NULL,
            mu0=list( trend=input$adjust.trend, #Trend adjustment?
                      S=season.adjust) #Seasonality? 0=no, 1 or 2 = # harmonics to include
          )) }else{
            mod1<- algo.glrnb(surv.ds1,control=list(
              range=c(n.weeks.train:input$week.test),
              c.ARL=input$set.thresh,
              M=-1, #How many time points back should we look? Negative 1: use all cases
              ret=c('value'),
              alpha=0,
              mu0=list( trend=input$adjust.trend, #Trend adjustment?
                        S=season.adjust) #Seasonality? 0=no, 1 or 2 = # harmonics to include
            ))
          }
        
        glr.vec<- c(rep(NA, times=(mod1$control$range[1]-1)),mod1$upperbound )
        m.vec<- c(rep(NA, times=(mod1$control$range[1]-1)),mod1$control$mu0 )
        
        alarm.vec<- c(rep(1, times=(mod1$control$range[1]-1)),mod1$alarm+2 )
        col.alarm.vec=c('gray', 'black','red')
        
        
        par(mfrow=c(2,1), mar=c(2,2,2,1))
        plot.obs <- mod1$disProgObj$observed[,1]
        plot.obs[(input$week.test+1):nrow(ds)] <- NA
        
        plot(ds[1:length(plot.obs),datevar],plot.obs, bty='l',type='p', ylab='Observed cases',col=col.alarm.vec[alarm.vec], pch=16)
        points(ds[1:length(m.vec),datevar],m.vec, type='l', col='black')
        title('Observed cases and mean')
        
        plot(ds[1:length(glr.vec),datevar],glr.vec, bty='l',type='p', ylab='GLR statistic',col=col.alarm.vec[alarm.vec], pch=16, xlim=c(min(ds[,datevar]),max(ds[,datevar]) ), 
             ylim=c(0, max(6,max(glr.vec, na.rm=T))))
        abline(h=mod1$control$c.ARL, lty=2)
        title('GLR statistic')
        
      },width = "auto", height = "auto")
    }
  )
}