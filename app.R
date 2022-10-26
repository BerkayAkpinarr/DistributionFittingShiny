library(shiny)
library(DT)
library(ggplot2)
library(cowplot)
library(scales)
library(dplyr)
library(tidyquant)
library(evir)
library(extRemes)
library(fitdistrplus)
library(EnvStats)
library(PearsonDS)
library(ggplot2)
library(reliaR)
library(goft)
library(gnFit)
library(lubridate)
library(readxl)
library(timetk)
library(Kendall)
library(trend)
##DATA SHOULD BE ORDERED AS "Year","Month","Day","Precipitation","Snow_Depth","Average_Temperature","Maximum_Temperature"


ui <- fluidPage(
  
  
  titlePanel(h1("Annual Based, Daily Precipitation Analysis with Fitting and Trend Analysis",align="center")),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose Your File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv",
                           ".xlsx"
                           )
                ),
      p("Please upload your data with the relevant headers:"),
      p("-Year, -Month, -Day, -Precipitation, -Snow_Depth, -Average_Temperature, -Maximum_Temperature", style = "font-family: 'times'; font-si16pt"),
      strong("Please be sure that your file has the mentioned headers"),
      p("What can be observed in Main Panels (Right Side):",style = "font-si16pt; color:blue"),
      p("Daily Max Values: This dataset created to obtain annual based daily maximum precipitation, occurence date and the snow depth at that time and maximum snow depth, occurence date and the precipitation at that time"),
      p("Fitted Values: This dataset gives fitted values for different distributions and return periods"),
      p("ECDF and Theoritical Plots: The plots show ECDF depictions of fitted distributions"),
      p("P values: The P-Values of Anderson-Darling and Cramer_von_Misses goodness of fit tests for distributions"),
      p("What can be observed in Lower Panels (Left/Below Side):",style = "font-si16pt; color:blue"),
      p("Trend Plots of Statistics of daily maximum precipitation and average temperature"),
      p("Quick and initial information of fitting with different plots"),
  
      fluidRow(
        column(h6(
          selectInput("se4", "Choices For Trend", choices = c("Nothing","Trend Plot","Trend Test"))
        ), width = 5, offset =
          0),
        column(h6(
          selectInput("se5", "Quick Fitting Information", choices = c("Nothing","GEV","GUM","LNORM","EXP"))
        ), width = 5, offset =
          0))
      ),
      
      #dateRangeInput("date_range","Date Range:",
      # start  = min(eddy_db_adjusted$TMSTAMP),
      # end    = max(eddy_db_adjusted$TMSTAMP),
      #textOutput("DateRange")),
     
    
        mainPanel(
          h3("Please do not forget to click -Nothing- when left hand side options used"),
          uiOutput("plots"),
          uiOutput("fitting_plots"),
           #h4(" ECDF and Theoritical Plots of Distributions"),
           #plotOutput("fitting_info_gev",width="900",height="600"),
          tabsetPanel(
            tabPanel(title= "Daily Max Values" ,DT::dataTableOutput("daily_max_dates")),
            tabPanel(title="Fitted Values",DT::dataTableOutput("fitting_values")),
            tabPanel(title="ECDF and Theoritical Plots",plotOutput("plot",width = 1000,height = 1500)),
            tabPanel(title="P Values",DT::dataTableOutput("p_values"))
          ),
      #DT::dataTableOutput("table")
    )
  )
)





server <- function(input, output) {
  
  output_data <- reactive({
    req(input$file1)
    df <- read_excel(input$file1$datapath)
    df$Date <- as.Date(with(df,paste(Year,Month,Day,sep="-")),"%Y-%m-%d")
    df <- df %>% relocate(Date,.after = Day)
    data_met <-df %>% timetk::summarise_by_time(
      .date_var = Date,
      .by = "year",
      Precipitation_Max = MAX(Precipitation),
      Precipitation_Sum = SUM(Precipitation),
      Temperature_Max = MAX(Maximum_Temperature),
      Temperature_Mean = AVERAGE(Average_Temperature),
      Max_Daily_Snow = MAX(Snow_Depth)
      
    )
    
    return(data_met)
    
  })
  
  
 
  output$fitting_values <- DT::renderDataTable({
    data_met <- output_data()
    p = c(0.995,0.99,0.98,0.95,0.9,0.8,0.5)
    gev.fit = fevd(data_met$Precipitation_Max,type ="GEV")
    loc <- as.numeric(gev.fit$results$par[1])
    scale <- as.numeric(gev.fit$results$par[2])
    shape <- as.numeric(gev.fit$results$par[3])
    gev.dist.results <- evd::qgev(p,loc,scale,shape)
    
    
    gum.fit = fevd(data_met$Precipitation_Max,type ="Gumbel")
    loc_gum <- as.numeric(gum.fit$results$par[1])
    scale_gum <- as.numeric(gum.fit$results$par[2])
    gumbel.dist.results <- evd::qgumbel(p,loc_gum,scale_gum)
   
    
    lnorm.fit <- fitdist(data_met$Precipitation_Max,"lnorm")
    lnorm_meanlog <- as.numeric(lnorm.fit$estimate[1])
    lnorm_sdlog <- as.numeric(lnorm.fit$estimate[2])
    lnormal.dist.results <- stats::qlnorm(p, lnorm_meanlog,lnorm_sdlog)
    
    
    exp.fit <- fitdist(data_met$Precipitation_Max,"exp")
    exp.dist.results <- stats::qexp(p,exp.fit$estimate)
    
    fitting_return_values <-rbind(gev.dist.results,gumbel.dist.results,lnormal.dist.results,exp.dist.results)
    colnames(fitting_return_values) <- c("200-yr","100-yr","50-yr","20-yr","10-yr","5-yr","2-yr")
    row.names(fitting_return_values) <- c("GEV","Gumbel","LogNormal","Exponential")
    fitting_return_values <- as.data.frame(fitting_return_values)
    
    return(fitting_return_values)
  })
  
  output$plot <- renderPlot({
    data_met <- output_data()
    p = c(0.995,0.99,0.98,0.95,0.9,0.8,0.5)
    gev.fit = fevd(data_met$Precipitation_Max,type ="GEV")
    loc <- as.numeric(gev.fit$results$par[1])
    scale <- as.numeric(gev.fit$results$par[2])
    shape <- as.numeric(gev.fit$results$par[3])
    gev.dist.results <- evd::qgev(p,loc,scale,shape)
    fitted <- evir::rgev(10000,shape,loc,scale)
    fitted_p <- evir::pgev(fitted,shape,loc,scale)
    fitted_gev = data.frame(fitted,fitted_p)
    
    gum.fit = fevd(data_met$Precipitation_Max,type ="Gumbel")
    loc_gum <- as.numeric(gum.fit$results$par[1])
    scale_gum <- as.numeric(gum.fit$results$par[2])
    gumbel.dist.results <- evd::qgumbel(p,loc_gum,scale_gum)
    fitted_gum <- evd::rgumbel(10000,loc_gum,scale_gum)
    fitted_gum_p <- evd::pgumbel(fitted_gum,loc_gum,scale_gum,lower.tail =T)
    fitted_gum_df = data.frame(fitted_gum,fitted_gum_p)
    
    lnorm.fit <- fitdist(data_met$Precipitation_Max,"lnorm")
    lnorm_meanlog <- as.numeric(lnorm.fit$estimate[1])
    lnorm_sdlog <- as.numeric(lnorm.fit$estimate[2])
    lnormal.dist.results <- stats::qlnorm(p, lnorm_meanlog,lnorm_sdlog)
    fitted_lnorm <- rlnorm(10000, lnorm_meanlog, lnorm_sdlog)
    fitted_lnorm_p <- plnorm(fitted_lnorm, lnorm_meanlog, lnorm_sdlog)
    fitted_lnorm_df = data.frame(fitted_lnorm,fitted_lnorm_p)
    
    exp.fit <- fitdist(data_met$Precipitation_Max,"exp")
    exp.dist.results <- stats::qexp(p,exp.fit$estimate)
    fitted_exp <- stats::rexp(10000, exp.fit$estimate)
    fitted_exp_p <- stats::pexp(fitted_exp, exp.fit$estimate)
    fitted_exp_df = data.frame(fitted_exp,fitted_exp_p)
    
    gev_plot <- ggplot(data_met,aes(Precipitation_Max)) + stat_ecdf(geom = "step",aes(colour = "ECDF"),show.legend = F)+
      geom_line(data = fitted_gev , aes(x=fitted, y =fitted_p), color = "blue") +labs(title = "GEV Distribution Fitting")+
      scale_x_continuous(breaks=pretty_breaks(n=8))+ 
      xlab("Maximum Precipitation (mm)") + ylab("Density") + scale_y_continuous(breaks=pretty_breaks(n=8))
      
    gum_plot <- ggplot(data_met,aes(Precipitation_Max)) + stat_ecdf(geom = "step",aes(colour = "ECDF"),show.legend = F)+
      geom_line(data = fitted_gum_df , aes(x=fitted_gum, y =fitted_gum_p), color = "blue") +
      scale_x_continuous(breaks=pretty_breaks(n=8))+
      labs(title = "Gumbel Distribution Fitting") +
      xlab("Maximum Precipitation (mm)") + ylab("Density") + scale_y_continuous(breaks=pretty_breaks(n=8))
    
    lnorm_plot <- ggplot(data_met,aes(Precipitation_Max)) + stat_ecdf(geom = "step",aes(colour = "ECDF"),show.legend = F)+
      geom_line(data = fitted_lnorm_df , aes(x=fitted_lnorm, y =fitted_lnorm_p), color = "blue") +
      scale_x_continuous(breaks=pretty_breaks(n=8))+
      labs(title = "Log-Normal Distribution Fitting") +
      xlab("Maximum Precipitation (mm)") + ylab("Density") + scale_y_continuous(breaks=pretty_breaks(n=8))
    
    exp_plot <- ggplot(data_met,aes(Precipitation_Max)) + stat_ecdf(geom = "step",aes(colour = "ECDF"),show.legend = F)+
      geom_line(data = fitted_exp_df , aes(x=fitted_exp, y =fitted_exp_p), color = "blue") +
      labs(title = "Exponential Distribution Fitting")+ scale_x_continuous(breaks=pretty_breaks(n=8))+
      xlab("Maximum Precipitation (mm)") + ylab("Density")+ scale_y_continuous(breaks=pretty_breaks(n=8))
  
    combined_plot <- plot_grid(gev_plot,gum_plot,lnorm_plot,exp_plot,
                             labels = "AUTO", 
                             ncol = 2,
                             align='hv',
                             label_fontfamily = "serif",
                             label_fontface = "plain",
                             label_colour = "blue",
                             label_x = 0, label_y = 0,
                             hjust = -0.5, vjust = -0.5)
    combined_plot
    #return(combined_plot)
    
  })
  output$daily_max_dates <- DT::renderDataTable({
    req(input$file1)
    df <- read_excel(input$file1$datapath)
    df$Date <- as.Date(with(df,paste(Year,Month,Day,sep="-")),"%Y-%m-%d")
    df <- df %>% relocate(Date,.after = Day)
    #df <- output_data()
    max_precip <-df %>% group_by(Year) %>% slice(which.max(Precipitation)) %>% dplyr::select(Date,Precipitation,Snow_Depth) %>% 
      rename(Max_Date_Precipitation = Date,  Maximum_Precipitation = Precipitation , Current_Snow_Depth = Snow_Depth)
    
    max_snow <-df %>%group_by(Year) %>% slice(which.max(Snow_Depth)) %>% dplyr::select(Date,Precipitation,Snow_Depth) %>% 
      rename(Max_Date_Snow = Date, Current_Precipitation = Precipitation  , Maximum_Snow_Depth = Snow_Depth)
    
    joined_max_values <- as.data.frame(full_join(max_precip,max_snow))
    return(joined_max_values)
  })
  
  
  output$p_values <- DT::renderDataTable({
    #req(input$file1)
    data_met <- output_data()
    lnorm_p_ad <-gofTest(data_met$Precipitation_Max, distribution = "lnorm",test = "ad")
    exp_p_ad <-gofTest(data_met$Precipitation_Max, distribution = "exp",test = "ad")
    lnorm_p_cvm <-gofTest(data_met$Precipitation_Max, distribution = "lnorm",test = "cvm")
    exp_p_cvm <-gofTest(data_met$Precipitation_Max, distribution = "exp",test = "cvm")
    
    gev.fit = fevd(data_met$Precipitation_Max,type ="GEV")
    gum.fit = fevd(data_met$Precipitation_Max,type ="Gumbel")
    gev_fit_p <- gnfit(data_met$Precipitation_Max,"gev",pr= gev.fit$results$par)
    gum_fit_p <- gnfit(data_met$Precipitation_Max,"gum",pr= gum.fit$results$par)
    
    p_values_df <- data.frame(Anderson_Darling = c(lnorm_p_ad$p.value,exp_p_ad$p.value,gev_fit_p$Apval,gum_fit_p$Apval),
                              Cramer_von_Misses = c(lnorm_p_cvm$p.value,exp_p_cvm$p.value,gev_fit_p$Wpval,gum_fit_p$Wpval))
    rownames(p_values_df) <- c("LogNormal","Exponential","GEV","Gumbel")
    
    return(p_values_df)  
    }) 
  output$plots <- renderUI({
   
    if(input$se4 == "Trend Plot"){
      h5(verbatimTextOutput("Plot"),
        plotOutput(
        "trend_plot", width = "1000px", height =
          "1000px"
      ), width = 1000)
    }else if(input$se4 == "Trend Test"){
      h5(dataTableOutput("trend_data"), width = 1000)
    }
  })
  
  output$trend_plot <- renderPlot({
    #req(input$file1)
    data_met <- output_data()
    data_met$Date <- year(data_met$Date)
    precipitation <- ggplot(data=data_met, aes(x=Date, y=Precipitation_Max)) +
      labs(title= "Precipitation vs Year")+theme_bw()+
      xlab("Year") + ylab("Precipitation (mm)")+
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=10,face="plain"))+
      scale_x_continuous(breaks = pretty_breaks(n=7))+
      scale_y_continuous(breaks = pretty_breaks(n=9))+
      geom_line(color="deepskyblue2")+
      geom_point(color="deepskyblue2")+
      geom_smooth(method=lm,color ="chocolate") #add linear trend line 
    
    temperature <-ggplot(data=data_met, aes(x=Date, y=Temperature_Mean, group=1)) +
      xlab("Year") + ylab("Temperature (???C)") + labs(title= "Average Temperature vs Year")+
      theme(axis.text=element_text(size=11),
            axis.title=element_text(size=11,face="plain"))+
      scale_x_continuous(breaks = pretty_breaks(n=7))+
      scale_x_continuous(breaks = pretty_breaks(n=9))+
      geom_line(color="deepskyblue2")+
      geom_point(color="deepskyblue2")+ theme_bw()+
      geom_smooth(method=lm,color ="chocolate") #add linear trend line 
    
    combined_plot <-  plot_grid(precipitation,temperature,labels = "AUTO", 
                                ncol = 1,
                                align='hv',
                                label_fontfamily = "serif",
                                label_fontface = "plain",
                                label_colour = "blue",
                                label_x = 0, label_y = 0,
                                hjust = -0.5, vjust = -0.5)
    return(combined_plot)
    
  })
  
  output$trend_data <- DT::renderDataTable({
    data_met <- output_data()
    kendall_precip <- MannKendall(data_met$Precipitation_Max)
    sens_slope_precip <- sens.slope(data_met$Precipitation_Max,conf.level = 0.95)
    
    kendall_temp <- MannKendall(data_met$Temperature_Mean)
    sens_slope_temp <- sens.slope(data_met$Temperature_Mean,conf.level = 0.95)
    
    trend_data <- data.frame(kendall_tau = c(kendall_precip$tau[1],kendall_temp$tau[1]),
                             kendall_p = c(kendall_precip$sl[1],kendall_temp$sl[1]),
                             sens_slope = c(sens_slope_precip$estimates[1],sens_slope_temp$estimates[1]),
                             sens_slope_p = c(sens_slope_precip$p.value,sens_slope_temp$p.value)
    )
    rownames(trend_data) <- c("Precipitation","Average Temperature")
    return(trend_data)
  })
 
  output$fitting_info_gev <- renderPlot({
    
    data_met <- output_data()
    gev.fit = fevd(data_met$Precipitation_Max,type ="GEV")
    plot(gev.fit)
  })
  
  output$fitting_info_gum <- renderPlot({
    data_met <- output_data()
    gum.fit = fevd(data_met$Precipitation_Max,type ="Gumbel")
    plot(gum.fit)
    
  })
  
  output$fitting_info_lnorm <- renderPlot({
    data_met <- output_data()
    lnorm.fit <- fitdist(data_met$Precipitation_Max,"lnorm")
    plot(lnorm.fit)
    
  })  
  
  output$fitting_info_exp <- renderPlot({
    data_met <- output_data()
    exp.fit <- fitdist(data_met$Precipitation_Max,"exp")
    plot(exp.fit)
  })  
  
  output$fitting_plots <- renderUI({
    
    if(input$se5 == "GEV"){
      h5(verbatimTextOutput("GEV Plot"),
         plotOutput(
           "fitting_info_gev", width = "1000px", height =
             "1000px"
         ), width = 1000)
    }else if(input$se5 == "GUM"){
      h5(verbatimTextOutput("GUM Plot"),
         plotOutput(
           "fitting_info_gum", width = "1000px", height =
             "1000px"
         ), width = 1000)
    }else if(input$se5 =="LNORM"){
      h5(verbatimTextOutput("LNORM Plot"),
         plotOutput(
           "fitting_info_lnorm", width = "1000px", height =
             "1000px"
         ), width = 1000)
    }else if(input$se5 =="EXP"){
      h5(verbatimTextOutput("EXP Plot"),
         plotOutput(
           "fitting_info_exp", width = "1000px", height =
             "1000px"
         ), width = 1000)
    }
    
  })
  
}    
 



shinyApp(ui = ui, server = server)



