# Graphs for presentation

# 1.Run 010_LifeExpSUI.R and 015_LifeDispSUI.R


# colorblind friendly pallete
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#F0E442", "#CC79A7")


# Graph 1 - Life Expectancy at age zero (male/female - faceted) #
# ------------------------------------------------------------- #

plotLE_zero <- LE %>% filter(Age==0) %>% ggplot() +
  # line plot
  geom_line(aes(x = Year, y = ex, color = cntry, alpha = highlight_flag))  +
  geom_point(aes(x = Year, y = ex, color = cntry, alpha = highlight_flag )) +
  scale_y_continuous(name = "Life expectancy at birth") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = c("black","grey20","grey35", "grey55","grey75","#FF6934"), name="") +
  scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  facet_wrap(~ sex) +
  theme_bw()
plotLE_zero <- plotLE_zero + theme(axis.text=element_text(size=12),
                                   axis.title=element_text(size=12,face="bold"))


# Graph 2 - Life Expectancy at age zero (males) #
# --------------------------------------------- #
plotLE_m_zero <- LE %>% filter(Age==0) %>% filter(sex=="male") %>% ggplot() +
  # line plot
  geom_line(aes(x = Year, y = ex, color = cntry, alpha = highlight_flag))  +
  geom_point(aes(x = Year, y = ex, color = cntry, alpha = highlight_flag )) +
  scale_y_continuous(name = "Life expectancy at birth") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = c("black","grey20","grey35", "grey55","grey75","#FF6934"), name="") +
  scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  theme_bw()
plotLE_m_zero <- plotLE_m_zero + theme(axis.text=element_text(size=12),
                                   axis.title=element_text(size=12,face="bold"))


# Graph 2 -  plot edagger trends #
# ------------------------------ #

edag_plot <- LSD %>% ggplot() +
  geom_line(aes(x = Year, y = ed_fem, color = cntry), linetype=1) + 
  geom_line(aes(x = Year, y = ed_mal, color = cntry), linetype=2) +
  scale_y_continuous(name = TeX('$e^\\dagger$')) +
  scale_colour_manual(values = cbbPalette, name="") +
  scale_alpha_discrete(range = c(0.35, 0.85), name="", guide=F) +  theme_bw()

edag_plot <- edag_plot + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))



# Graph 3 - Plot the (female - male) gap in life expectancy at birth #
# ------------------------------------------------------------------ #



plotgap <- LE_GAP %>% ggplot() +
  # line plot
  geom_line(aes(x = Year, y = gap, color=cntry))  +
  geom_point(aes(x = Year, y = gap, color=cntry)) +
  scale_y_continuous(name = "Female-Male Gap in LE at birth") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = cbbPalette, name="", guide=F) +
  scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  theme_bw()

plotgap <- plotgap +  scale_shape_discrete(guide=FALSE) + theme(axis.text=element_text(size=12),
                                                                        axis.title=element_text(size=12,face="bold"))

# Graph 4 - Plot the (female - male) gap in life expectancy at birth #
# ------------------------------------------------------------------ #
# calculate the gap in LE
LSD <- LSD %>% mutate(gap = ed_mal - ed_fem)

# Life Expectancy at age 65
plotgap_LSD <- LSD %>% ggplot() +
  # line plot
  geom_line(aes(x = Year, y = gap, color=cntry))  +
  geom_point(aes(x = Year, y = gap, color=cntry)) +
  scale_y_continuous(name = "Male-Female Gap in Life Span Disparity") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = cbbPalette, name="", guide=F) +
  scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  theme_bw()

plotgap_LSD <- plotgap_LSD +  scale_shape_discrete(guide=FALSE) + theme(axis.text=element_text(size=12),
                                                                                axis.title=element_text(size=12,face="bold"))

multiplot(plotgap, plotgap_LSD)

# Graph 5 - Plot the (female - male) gap in life expectancy at birth / Age 65 - highlighting Switzerland
# ------------------------------------------------------------------------------------------------------ #

plotgap_SUI <- LE_GAP %>% ggplot() +
  # line plot
  geom_line(aes(x = Year, y = gap, color=cntry, alpha = highlight_flag))  +
  geom_point(aes(x = Year, y = gap, color=cntry, alpha = highlight_flag)) +
  scale_y_continuous(name = "Female-Male Gap in LE at birth") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = c("black","grey20","grey35", "grey55","grey75","#FF6934"), name="", guide=F) +
  scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  theme_bw()

plotgap_SUI <- plotgap_SUI +  scale_shape_discrete(guide=FALSE) + theme(axis.text=element_text(size=12),
                                                                axis.title=element_text(size=12,face="bold"))

# Graph 6 - Plot the (female - male) gap in life span disparity (e dagger) - highlighting Switzerland #
# --------------------------------------------------------------------------------------------------- #
# calculate the gap in LE
LSD <- LSD %>% mutate(gap = ed_mal - ed_fem)

# Life Expectancy at age 65
plotgap_LSD_SUI <- LSD %>% ggplot() +
  # line plot
  geom_line(aes(x = Year, y = gap, color=cntry, alpha = highlight_flag))  +
  geom_point(aes(x = Year, y = gap, color=cntry, alpha = highlight_flag)) +
  scale_y_continuous(name = "Male-Female Gap in Life Span Disparity") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = c("black","grey20","grey35", "grey55","grey75","#FF6934"), name="", guide=F) +
  scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  theme_bw()

plotgap_LSD_SUI <- plotgap_LSD_SUI +  scale_shape_discrete(guide=FALSE) + theme(axis.text=element_text(size=12),
                                                                        axis.title=element_text(size=12,face="bold"))

# 6 a-b Plot both together

# run multiplot code

# plot again
multiplot(plotgap_SUI,plotgap_LSD_SUI)



# Graph 7-8 - Plot the (female - male) gap in life expectancy at birth / e-dagger - highlighting Denmark #
# ------------------------------------------------------------------------------------------------------ #

LE_GAP <- LE_GAP %>% mutate(highlight_flag_2 = ifelse(cntry=="DEN",T,F))

plotgap_DEN <- LE_GAP %>% ggplot() +
  # line plot
  geom_line(aes(x = Year, y = gap, color=cntry, alpha = highlight_flag_2))  +
  geom_point(aes(x = Year, y = gap, color=cntry, alpha = highlight_flag_2)) +
  scale_y_continuous(name = "Female-Male Gap in LE at birth") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = c("#009E73","grey20","grey35", "grey55","grey75","black"), name="", guide=F) +
  scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  theme_bw()

plotgap_DEN <- plotgap_DEN +  scale_shape_discrete(guide=FALSE) + theme(axis.text=element_text(size=12),
                                                                        axis.title=element_text(size=12,face="bold"))

# calculate the gap in LE
LSD <- LSD %>% mutate(highlight_flag_2 = ifelse(cntry=="DEN",T,F))

# Life Expectancy at age 65
plotgap_LSD_DEN <- LSD %>% ggplot() +
  # line plot
  geom_line(aes(x = Year, y = gap, color=cntry, alpha = highlight_flag_2))  +
  geom_point(aes(x = Year, y = gap, color=cntry, alpha = highlight_flag_2)) +
  scale_y_continuous(name = "Male-Female Gap in Life Span Disparity") +
  scale_x_continuous(name = " ") +
  scale_colour_manual(values = c("#009E73","grey20","grey35", "grey55","grey75","black"), name="", guide=F) +
  scale_alpha_discrete(range = c(0.25, 0.85), name="", guide=F) +
  theme_bw()

plotgap_LSD_DEN <- plotgap_LSD_DEN +  scale_shape_discrete(guide=FALSE) + theme(axis.text=element_text(size=12),
                                                                                axis.title=element_text(size=12,face="bold"))

# plot again
multiplot(plotgap_DEN,plotgap_LSD_DEN)

# "grey20","grey35", "grey55","grey75",
# "#26B7FF","#26FF57", "#FFD846","#FF6600","#0D3BB2","#FF6934"