### Thanks to metafor-project.org 
### https://www.metafor-project.org/doku.php/plots:forest_plot_revman

### You will want to use png or pdf or tif to save this res 350, width=3196, height=1648

# Ensures the package "pacman" is installed
if (!require("pacman")) install.packages("pacman")

### You will also need cmdstanr and rstan installed
### If you are on R 4.2+ you'll need R tools to install them as well as this link
### https://blog.mc-stan.org/2022/04/26/stan-r-4-2-on-windows/

### Import Dataset
pacman::p_load(readxl,
               metafor,
               brms,
               dplyr,
               tidybayes,
               ggplot2,
               bayesmeta)

#Create the text for tau
tau_text <- function(model) {
  
  posterior = 
    model |> 
    tidybayes::tidy_draws() |> 
    ggdist::median_hdi(sd_author__Intercept)
  
  bquote(paste("Heterogeneity: ",
               
               tau,
               .(" = "),
               .(formatC(posterior$sd_author__Intercept, digits=2, format="f")),
               .(" [95% CrI: "),
               .(formatC(posterior$.lower, digits=2, format="f")),
               ", ",
               .(formatC(posterior$.upper, digits=2, format="f")),
               "]"
  )
  )
}


#Read the data
dat <- read_excel("data/data.xlsx")

### turn the risk of bias into levels
dat[7:12] <- lapply(dat[7:12], factor, levels=c("+", "-", "?"))

### calculate log risk ratios and corresponding sampling variances (and use
### the 'slab' argument to store study labels as part of the data frame)

dat <- escalc(measure="RR", ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat,
              slab=paste(author, year), drop00=TRUE)

#Add sei from vi
dat$sei = sqrt(dat$vi)

# Bayes

## Formula
mf = 
  # https://bookdown.org/content/4857/horoscopes-insights.html#consider-using-the-0-intercept-syntax
  formula(yi | se(sei) ~ 0 + Intercept + (1 | author))

## Priors

# Tau is based on Turner Et Al https://pubmed.ncbi.nlm.nih.gov/25475839/
informative = bayesmeta::TurnerEtAlPrior("all-cause mortality",
                                         "pharmacological",
                                         "placebo / control")


logmean = informative$parameters["tau", "mu"]
logsd = informative$parameters["tau", "sigma"]

# logmean
# -1.975

# logsd
# 0.67

#Prepare the priors; here we use a weakly informative prior centred at RR 1 with 95% CI between 0.25 and 4
#This is because the majority of pharmaceutical interventions which work will be within this range
priors_ma = 
  brms::prior(normal(0, 0.71), class = "b", coef = "Intercept") +
  brms::prior(lognormal(-1.975, 0.67), class = "sd") 

## Models

ma_bayes = 
  brms::brm(
    data = dat,
    family = gaussian,
    
    formula = mf,
    prior = priors_ma,
    sample_prior = TRUE,
    
    control = list(adapt_delta = .95),
    backend = "cmdstanr", # faster
    cores = parallel::detectCores(),
    chains = 4,
    warmup = 5000, 
    iter = 10000, 
    seed = 123,
    file = "output/models/ma_bayes.Rds",
    file_refit = "on_change"
  )

ma_bayes

plot(ma_bayes)

pp_check(ma_bayes)

#This function makes the master figure which works best using PNG output as commented below and not X11/RStudio window
forest_plot = function(){
  
  ### estimated average risk ratio
  
  pred = ma_bayes |> 
    brms::fixef(summary = F) |>
    ggdist::median_hdi() 
  
  predt = pred |> 
    mutate(across(y:ymax, ~exp(.)))
  
  ### need the rounded estimate and CI bounds further below
  predt <- formatC(c(predt$y, predt$ymin, predt$ymax), format="f", digits=2)
  
  ### total number of studies
  k <- nrow(dat)
  
  ### set na.action to "na.pass" (instead of the default, which is "na.omit"),
  ### so that even a study with missing log odds/risk ratio will be shown in the
  ### forest plot
  
  options(na.action = "na.pass")
  
  ### adjust the margins
  par(mar=c(12,0,1.3,1.3), mgp=c(3,0.2,0), tcl=-0.2)
  
  ### forest plot with extra annotations
  
  sav <- with(dat,
              forest(yi, vi,
                     atransf=exp,
                     at=log(c(0.01,.10, 0.25,0.5, 1, 2,4,10,100)),
                     xlim=c(-30,11),
                     ylim = c(0, 10),
                     xlab="", efac=c(0,4), textpos=c(-30,-4.7), lty=c(1,1,0), 
                     refline=NA,
                     ilab=cbind(ai, n1i, ci, n2i),
                     ilab.xpos=c(-20.6,-18.6,-16.1,-14.1), ilab.pos=2,
                     cex=0.78, header=c("Study","Median, 95% CI"), mlab="")
  )
  
### add horizontal line at the top
  segments(sav$xlim[1]+0.5, k+1, sav$xlim[2], k+1, lwd=0.8)
  
  ### add vertical reference line at 0
  segments(0, -2, 0, k+1, lwd=0.8)
  
  ### now we add a bunch of text; since some of the text falls outside of the
  ### plot region, we set xpd=NA so nothing gets clipped
  par(xpd=NA)
  
  ### adjust cex as used in the forest plot and use a bold font
  par(cex=sav$cex, font=2)
  
  #Add headers
  text(sav$ilab.xpos, k+2, pos=2, c("Events","Total","Events","Total"))
  text(c(mean(sav$ilab.xpos[1:2]),mean(sav$ilab.xpos[3:4])), k+3, pos=2,
       c("Corticosteroids","Control"))
  #text(sav$textpos[2], k+3, "Risk Ratio", pos=2)
  text(0, k+2, "Risk ratio")
  text(sav$xlim[2]-0.6, k+3, "Risk of Bias", pos=2)
  text(c(sav$xlim[1],sav$ilab.xpos[c(2,4,5)]), 0, pos=c(4,2,2,2,2),
       c("Total (95% CrI)", sum(dat$n1i), sum(dat$n2i)))
  ### add total events
  text(sav$ilab.xpos[c(1,3)], 0, c(sum(dat$ai),sum(dat$ci)), pos=2)
  
  #Add header for Risk of bias legend
  text(sav$xlim[1], -3, pos=4, "Risk of bias legend")
  
  ### use a non-bold font for the rest of the text
  par(cex=sav$cex, font=1)
  
  ### add 'Favours '/'Favours proph' text below the x-axis scale
  text(log(c(.01, 100)), -2, c("Favors Corticosteroids","Favors Control"),
       pos=c(4,2), offset=-0.5)
  
  #Add the heterogeneity text using the tau_text function
  text(sav$xlim[1], -1,
       tau_text(ma_bayes),
       pos=4,
       cex=1)
  #Add a diamond for the overall analysis based on the predicted values from the bayesian CrI
  metafor::addpoly(x = pred$y,
                   ci.lb = pred$ymin,
                   ci.ub = pred$ymax,
                   
                  # mlab = "Total [95% CrI]", # label
                   row=0, # location
                   efac = 3 # polygon size
  )
  
  ### first hide the non-bold summary estimate text and then add it back in bold font
  rect(sav$textpos[2] - 4, -0.5,
       sav$textpos[2] , 0.5,
       col="white", border=NA)
  

  #go bold!
  par(cex=sav$cex, font=2)
  # Now replace the text
  text(sav$textpos[2], 0, paste0(predt[1], " [", predt[2],
                                 ",  ", predt[3], "]"),
       pos=2, bold = 2)
  # Stop bold
  par(cex=sav$cex, font=1)
  
  ### add risk of bias points and symbols
  cols <- c("#00cc00", "#cc0000", "#eeee00")
  syms <- levels(dat$rb.a)
  pos  <- seq(sav$xlim[2]-5.5,sav$xlim[2]-0.5,length=6)
  for (i in 1:6) {
    points(rep(pos[i],k), k:1, pch=19, col=cols[dat[[6+i]]], cex=2.2)
    text(pos[i], k:1, syms[dat[[6+i]]], font=2)
  }
  text(pos, k+2, c("A","B","C","D","E","F"), font=2)
  
  ### add risk of bias legend below the heading above
  text(sav$xlim[1], -4:-9, pos=4, c(
    "(A) Bias arising from the randomization process",
    "(B) Bias due to deviations from the intended intervention",
    "(C) Bias due to missing outcome data",
    "(D) Bias in measurement of the outcome",
    "(E) Bias in selection of the reported result",
    "(F) Overall"))
  
}

#png(filename="C:/output/forest.png", res=350, width=3196, height=1648)
forest_plot()
#dev.off()
#x11()

###Probability plot -

#Here we set prob2 based on being below the RR corresponding to 4% ARR from weighted control event rate of 18%
#Prob 1 implies below 0 - e.g., ARR>=0%
probs = 
  ma_bayes |> 
  tidy_draws() |> 
  summarise(prob2 = 100*mean(b_Intercept < log(14/18)),
            prob1 = 100*mean(b_Intercept < log(1)))

#Create the textbox
text1<-paste0("Probability ARR >0% (dark + light purple): ",
              round(probs$prob1,1) ,
              "%\nProbability ARR >=4% (light purple): ",
              round(probs$prob2,1) |> paste0("%"))

#Create a table we can reference in the graph
overall = fixef(ma_bayes) |> data.frame()


#png(filename="C:/output/probability.png", res=350, width=3196, height=1648)
ggplot(data.frame(x = c(log(0.25), log(1.5))), aes(x = x)) +
  geom_vline(xintercept = 0, linetype="dashed",color="red")+  
  theme_bw()+
  stat_function(fun = dnorm, args=c(overall$Estimate, overall$Est.Error),
                size=0.5)  +
  geom_area(stat ="function", fun=dnorm, fill="darkblue", alpha=0.2, xlim=c(log(0.25),log(14/18)),
            args=c(overall$Estimate, overall$Est.Error))+
  geom_area(stat ="function", fun=dnorm, fill="darkblue", alpha=0.6, xlim=c(log(14/18),0),
            args=c(overall$Estimate, overall$Est.Error))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major=element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(angle = 180)) +
  scale_y_continuous(position="right",
                     expand = c(0,0))+
  labs(
    x="Risk Ratio (log scale)",
    y="Density\n"
  )+
  geom_label(aes(x=-1,y=2.25,label=text1),size=4) +
  scale_x_continuous(labels = c(0.25,seq(0.5, 1.25, 0.25)), 
                     breaks = log(c(0.25, seq(0.5, 1.25, 0.25)))) +
  coord_cartesian(x = log(c(0.25, 1.3)),
                  y = c(0, 2.5))

#dev.off()
#x11()
