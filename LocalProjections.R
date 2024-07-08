#####################
#Kevin Corfield
#Universidad de Buenos Aires 
#Local Projections
#####################
library(lmtest)
library(sandwich)
library(ivreg)
library(dynlm)
library(haven)
library(ggplot2)
data_monetary_shock_quarterly <- read_stata("C:/Users/Kevin/Desktop/LPIV/data_monetary_shock_quarterly.dta")
View(data_monetary_shock_quarterly)

# Keep only nonmissing observations in resid_full i.e. 1969q1 - 2007q4
df <- subset(data_monetary_shock_quarterly, !is.na(resid_full))
df <- df[,-1]

#Choose impulse response horizon
hmax=16

# Generate LHS variables for the LPs (Cumulative)
for (i in 1:hmax+1) {
  nm <- paste("ur", i-1, sep = "")
  df[[nm]] <- c(df$UNRATE[i:nrow(df)], rep(NA, i-1))
}


#########################################
# Convert to Tseries object
#########################################
library(tseries)
df_ts <- ts(df, start = c(1969,1), frequency = 4)
View(df)

#################################################
#Local projections by OLS Using Funds Rates
#################################################

bh <- matrix(nrow = hmax, ncol = 3)

for (h in 1:hmax) {
      
      reg <- dynlm(df_ts[,ncol(df_ts)-hmax+h] ~ DFF + L(DFF,1) + L(DFF,2) + L(DFF,3) + L(DFF,4) 
                    + L(UNRATE,1) + L(UNRATE,2) + L(UNRATE,3) + L(UNRATE,4), data = df_ts)  
      
      coefs=coeftest(reg, vcov=NeweyWest(reg, lag = h, prewhite = FALSE, adjust = T))
        
        #Store coefficient of DFF
        bh[h,1] = reg$coefficients[[2]]
        
        # compute standard error
        #sd <- sqrt(diag(NW_VCOV))[[2]]
        #se <- sd/sqrt(length(UNRATE)-h)
        
        #Store Upper and Lower bands 90% confidence Bands
        bh[h,2] <- bh[h,1] + 1.645 * coefs[2,2]
        bh[h,3] <- bh[h,1] - 1.645 * coefs[2,2]
}

#Plot results
x=1:hmax

bh <-data.frame(bh)
(p <- ggplot(bh, aes(x, bh[,1]))+
    geom_line(data=bh, size=2)+
    ggtitle("Responses of the unemployment rate to monetary shock") +
    scale_color_manual(values="purple")+
    geom_ribbon(data=bh,aes(ymin=bh[,2],ymax=bh[,3]),alpha=0.2, fill = "purple")+
    theme_light() +
    xlab("Quarter")+
    ylab("Percent"))


#######################################################################################
#Local projections by IV
#######################################################################################
bh_iv <- matrix(nrow = hmax, ncol = 3)

for (h in 1:hmax) {
  
  reg_iv <- ivreg(df_ts[,ncol(df_ts)-hmax+h] ~ DFF + lag(DFF,1) + lag(DFF,2) + lag(DFF,3) + lag(DFF,4) 
                  + lag(UNRATE,1) + lag(UNRATE,2) + lag(UNRATE,3) + lag(UNRATE,4) | resid_full +  
                    lag(resid_full,1) + lag(resid_full,2) + lag(resid_full,3) + lag(resid_full,4)  + 
                    lag(UNRATE,2) + lag(UNRATE,3) + lag(UNRATE,4), 
                    data = df_ts)  
  
  coefs=coeftest(reg_iv, vcov=NeweyWest(reg_iv, lag = h, prewhite = FALSE, adjust = T))
  
  #Store coefficient of DFF
  bh_iv[h,1] = reg_iv$coefficients[[2]]
  
  # compute standard error
  #sd <- sqrt(diag(NW_VCOV))[[2]]
  #se <- sd/sqrt(length(UNRATE)-h)
  
  #Store Upper and Lower bands
  bh_iv[h,2] <- bh_iv[h,1] + 1.645 * coefs[2,2]
  bh_iv[h,3] <- bh_iv[h,1] - 1.645 * coefs[2,2]
}

#######################################################################################
# Plot
#######################################################################################
bh_iv <-data.frame(bh_iv)
# Plot bh_iv
plot_bh_iv <- ggplot(bh_iv, aes(x, bh_iv[,1])) +
  geom_line(size = 2, aes(color = "purple")) +
  geom_ribbon(aes(ymin = bh_iv[,2], ymax = bh_iv[,3]), alpha = 0.2, fill = "purple") +
  ggtitle("Responses of the unemployment rate to monetary shock") +
  xlab("Quarter") +
  ylab("Percent") +
  theme_light() +
  geom_hline(yintercept = 0, linetype = "dashed")  

# Plot bh
plot_bh <- ggplot(bh, aes(x, bh[,1])) +
  geom_line(size = 2, aes(color = "blue"), linetype = "dashed") +  
  geom_ribbon(aes(ymin = bh[,2], ymax = bh[,3]), alpha = 0.2, fill = "blue") +
  xlab("Quarter") +
  ylab("Percent") +
  theme_light() +
  geom_hline(yintercept = 0, linetype = "dashed")  

# Combinar los gráficos
combined_plot <- plot_bh_iv +
  geom_line(data = bh, aes(x, bh[,1], color = "blue")) +
  geom_ribbon(data = bh, aes(x, ymin = bh[,2], ymax = bh[,3]), alpha = 0.2, fill = "blue") +
  scale_color_manual(name = "", values = c("purple" = "purple", "blue" = "blue"),
                     labels = c("OLS", "IV")) +
  guides(color = guide_legend(title = "Estimation")) +  
  geom_hline(yintercept = 0, linetype = "solid") +  
  theme(legend.position = c(0.1, 0.9), legend.justification = c(0, 1),
        legend.background = element_rect(fill = "transparent"),  
        legend.key = element_rect(color = NA, fill = "transparent"))  


# Mostrar el gráfico combinado
print(combined_plot)
