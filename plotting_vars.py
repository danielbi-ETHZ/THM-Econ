## For plotting
## refer to https://www.rapidtables.com/web/color/brown-color.html
# blanchedalmond = (255/255.,235/255.,205/255.)
navajowhite = (255/255.,222/255.,173/255.)
burlybrown = (222/255.,184/255.,135/255.)
peru = (205/255.,133/255.,64/255.)
saddlebrown = (139/255.,69/255.,19/255.)
maroon = (128/255.,30/255.,15/255.)
colors = [navajowhite, burlybrown, peru, saddlebrown]

### colors from https://github.com/OrdnanceSurvey/GeoDataViz-Toolkit/blob/master/Colours/GDV%20colour%20palettes%200.7.pdf
Red = '#FF1F5B'
Purple = '#AF58BA'
Blue = '#009ADE'
Green = '#00CD6C'
Yellow = '#FFC6aE'

### For contour plot
good_green = (38/256., 153/256., 0/256.)
ok_yellow = (179/256., 179/256., 0/256.)
ok_red = '#852810' # (87/256., 18/256., 0/256.)
# bad_red = (179/256., 36/256., 0/256.)
# worse_red = (204/256., 41/256., 0/256.)

### Red gradations for fig7a.py
color200 = (87/256., 18/256., 0/256.)# ok_red
color500 = ok_red #'#852810'
color800 = '#ad4328'
color1500 = '#d16145'
color2500 = '#ed876d'

### Green for Figure 4(c) (Lstar_vs_k_fig2b.py)
color200_green = '#06592A'
color500_green = '#228B3B'
color800_green = '#40AD5A'
color1500_green = '#6CBA7D'
color2500_green = '#9CCEA7'


normal_width = 1
base_case_width = 2
normal_font_size = 12 ## 13 ## for (a), (b) labels
label_font_size = 14 ## for xlabel, ylabel
tick_font_size =12 # 12 ## for 0, 10, 20...
marker_size = 9

#### Define the d^*_{LCOH,min}, read from plot of depth vs LCOH
dstar_LCOH_min_alpha1 = [-460,-273,-189] # depths of min LCOH for b = [10,20,40], alpha = 1.0
dstar_LCOH_min_alpha08 = [-742,-392,-273] # depths of min LCOH for b = [10,20,40], alpha = 0.8
#### Define the d^*_{kmin}, read from plot of depth versus kmin
dstar_kmin_alpha1 = [-496,-465,-453] # the depth of the minimum kmin, for alpha = 1.0, b = [10,20,40]
dstar_kmin_alpha08 = [-692,-563,-528] # the depth of the minimum kmin, for alpha = 0.8, b = [10,20,40]

FIGURE_DIRECTORY = '.'

