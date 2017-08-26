# Binodal curve at different temperature
rm(list = ls())
options(digits = 10)
source('vomodel.R')
source('para.voorn.chi.R')
library(yxplot)
library(ggplot2)
library(dplyr)
library(rootSolve)
library(pracma)
library(nleqslv)
DEBUG <- T
SAVE <- F

fitting.para$sampling.end <- 1e-3
function(){
}
test.chi.chipw <- function(DEBUG = T, SAVE = F) {
        
    fitting.para$binodal.guess <- c(0.05, 0.03)
    # fitting.para$sampling.start <- 0.0007875
    
    ds <- NULL
    # Chi protein and water
    for (chipw in c( -0.01, 0, 0.01)) {
        system.properties$Chi <- matrix(c(.0,.0,0,0, chipw,
                                        .0,.0,0,0,0,
                                        0,0,0,0,0,
                                        0,0,0,0,0,
                                        chipw,0,0,0,0), 5,5)
        
        p20 <- get.binodal.curve(80, system.properties$Chi, system.properties, fitting.para) 
        if (is.null(p20)) return(NULL)
        p20 <- p20 %>% mutate(Chi = chipw)
        ds <- rbind(ds, p20)
    }
    
    # critical points
    print(unique(ds$critic.salt))
    # Plots
    g2 <- ggplot(ds, aes(y = round(phi.salt*100, digits = 4))) +
        geom_point(aes(x = phi.polymer*100, col = as.factor(Chi)), lwd = 2) +
        # geom_line(aes(x = phi.polymer*100, group = pairing, col = as.factor(Chi))) +
        labs(x = 'Polymer [%]',
             y = 'Salt [%]',
             col = 'Chi pw')
    g2 <- theme.title.text.1(g2)
    
    print(g2)
    if(SAVE) {
        ggsave('test.chi.chipw.png', g2, width = 6, height = 3)
    } 
}
test.chi.chipw(DEBUG, SAVE)
test.chi.chipq <- function(DEBUG = T, SAVE = F) {
    
fitting.para$binodal.guess <- c(0.03, 0.05)
fitting.para$critical.point.guess <- c(0.01, 0.001)
# fitting.para$sampling.start <- 0.0007875

ds <- NULL
# Chi protein and water
for (chipq in c(-0.05, 0, 0.05)) {
system.properties$Chi <- matrix(c(0,chipq,0,0, 0,
                                chipq,.0,0,0,0,
                                0,0,0,0,0,
                                0,0,0,0,0,
                                0,0,0,0,0), 5,5)
p20 <- get.binodal.curve(80, system.properties$Chi, system.properties, fitting.para) %>% 
    mutate(Chi = chipq)
ds <- rbind(ds, p20)
}

# critical points
print(unique(ds$critic.salt))
# Plots
g2 <- ggplot(ds, aes(y = round(phi.salt*100, digits = 4))) +
    geom_point(aes(x = phi.polymer*100, col = as.factor(Chi)), lwd = 2) +
    # geom_line(aes(x = phi.polymer*100, group = pairing, col = as.factor(Chi))) +
    labs(x = 'Polymer [%]',
         y = 'Salt [%]',
         col = 'Chi pq')
g2 <- theme.title.text.1(g2)


print(g2)
if(SAVE) {
    ggsave('test.chi.chipq.png', g2, width = 6, height = 3)
} 
}
test.chi.chipq(DEBUG, SAVE)
test.chi.chips <- function(DEBUG = T, SAVE = F) {
    
    fitting.para$binodal.guess <- c(0.03, 0.05)
    fitting.para$critical.point.guess <- c(0.01, 0.001)
    # fitting.para$sampling.start <- 0.0007875
    
    ds <- NULL
    # Chi protein and water
    for (chips in c(-0.05, 0, 0.05)) {
    system.properties$Chi <- matrix(c(
        0,    0,    chips,       0,   0,
        0,   .0,        0,   chips,   0,
        chips, 0,       0,       0,   0,
        0,    chips,    0,       0,   0,
        0,     0,       0,       0,   0), 5,5)
    p20 <- get.binodal.curve(80, system.properties$Chi, system.properties, fitting.para) %>% 
        mutate(Chi = chips)
    ds <- rbind(ds, p20)
    }
    
    # critical points
    print(unique(ds$critic.salt))
    # Plots
    g2 <- ggplot(ds, aes(y = round(phi.salt*100, digits = 4))) +
        geom_point(aes(x = phi.polymer*100, col = as.factor(Chi)), lwd = 2) +
        # geom_line(aes(x = phi.polymer*100, group = pairing, col = as.factor(Chi))) +
        labs(x = 'Polymer [%]',
             y = 'Salt [%]',
             col = 'Chi ps')
    g2 <- theme.title.text.1(g2)
    
    
    print(g2)
    if(SAVE) {
        ggsave('test.chi.chips.png', g2, width = 6, height = 3)
    } 
}
test.chi.chips(DEBUG, SAVE)

function(){
    
test.chipw.temp <- function(DEBUG = T, SAVE = F) {
    
fitting.para$binodal.guess <- c(0.05, 0.1)
# fitting.para$sampling.start <- 0.0007875

ds <- NULL
# Chi protein and water
chipw.0 <- -0.01
for (temp in c(20, 30)) {
    
    chipw <- chipw.0 * (20 + 273) / (temp+273)
    
    system.properties$Chi <- matrix(c(.0,.0,0,0, chipw,
                                    .0,.0,0,0,0,
                                    0,0,0,0,0,
                                    0,0,0,0,0,
                                    chipw,0,0,0,0), 5,5)
    p20 <- get.binodal.curve(temp, system.properties$Chi, system.properties, fitting.para) %>% 
        mutate(Chi = chipw)
    ds <- rbind(ds, p20)
    }
    
    # critical points
    print(unique(ds$critic.salt))
    # Plots
    g2 <- ggplot(ds, aes(y = phi.salt*100)) +
        geom_point(aes(x = phi.polymer*100, col = as.factor(tempC)), lwd = 2) +
        # geom_line(aes(x = phi.polymer*100, group = pairing, col = as.factor(Chi))) +
        labs(x = 'Polymer [%]',
             y = 'Salt [%]',
             col = 'TempC')
    g2 <- theme.title.text.1(g2)
    
    
    print(g2)
    if(SAVE) {
        ggsave('test.binodal.curve.voorn.png', g2, width = 6, height = 3)
    } 
}
test.chipw.temp()
test.chi.chipq.temp <- function(DEBUG = T, SAVE = F) {
    
fitting.para$binodal.guess <- c(0.03, 0.05)
fitting.para$critical.point.guess <- c(0.01, 0.0001)
# fitting.para$sampling.start <- 0.0007875

ds <- NULL
# Chi protein and water
chipq.0 <- -0.05
for (temp in c(20, 30)) {
    chipq <- chipq.0 * (80 + 273) / (temp+273)
    # Chi
    print(c('Chi', chipq))
        
    system.properties$Chi <- matrix(c(0,chipq,0,0, 0,
                                    chipq,.0,0,0,0,
                                    0,0,0,0,0,
                                    0,0,0,0,0,
                                    0,0,0,0,0), 5,5)
    p20 <- get.binodal.curve(temp, system.properties$Chi, system.properties, fitting.para) %>% 
        mutate(Chi = chipq)
    ds <- rbind(ds, p20)

}

# critical points
print(unique(ds$critic.salt))

# Plots
g2 <- ggplot(ds, aes(y = round(phi.salt*100, digits = 4))) +
    geom_point(aes(x = phi.polymer*100, col = as.factor(tempC)), lwd = 2) +
    # geom_line(aes(x = phi.polymer*100, group = pairing, col = as.factor(Chi))) +
    labs(x = 'Polymer [%]',
         y = 'Salt [%]',
         col = 'TempC')
g2 <- theme.title.text.1(g2)


print(g2)
if(SAVE) {
    ggsave('test.chi.chipq.png', g2, width = 6, height = 3)
} 
}
test.chi.chipq.temp()
}
