import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
import matplotlib as mp



plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
mp.rcParams["font.size"] = "12.5"

fig, axs = plt.subplots(1, 2, figsize = (14.5, 5))


source = Rectangle((1 - 0.25, 1), 0.5, 8, facecolor = 'k')
etalon = Rectangle((6 - 0.25, 1), 0.5, 8, facecolor = 'k')
lens = Rectangle((11 - 0.25, 1), 0.5, 8, facecolor = 'k')
sensor = Rectangle((19 - 0.25, 1), 0.5, 8, facecolor = 'k')


lr1 = axs[0].annotate("", xy=(10.5, 8), xytext=(1.5, 2),arrowprops=dict(arrowstyle="->", color = 'k', alpha = 0.5))
lr2 = axs[0].annotate("", xy=(10.5, 5), xytext=(1.5, 5),arrowprops=dict(arrowstyle="->", color = 'k', alpha = 0.5))
lr3 = axs[0].annotate("", xy=(10.5, 2), xytext=(1.5, 8),arrowprops=dict(arrowstyle="->", color = 'k', alpha = 0.5))


flr1 = axs[0].annotate("", xy=(18.5, 8), xytext=(11.5, 8),arrowprops=dict(arrowstyle="->", color = 'k', alpha = 0.5))
flr2 = axs[0].annotate("", xy=(18.5, 5), xytext=(11.5, 5),arrowprops=dict(arrowstyle="->", color = 'k', alpha = 0.5))
flr3 = axs[0].annotate("", xy=(18.5, 2), xytext=(11.5, 2),arrowprops=dict(arrowstyle="->", color = 'k', alpha = 0.5))


width = 0.5
axs[0].fill_between([1.5, 10.5], [2 - width, 8 - width], [2 + width, 8 + width], facecolor = 'deeppink', alpha = 0.3)
axs[0].fill_between([1.5, 10.5], [8 - width, 2 - width], [8 + width, 2 + width], facecolor = 'darkorange', alpha = 0.3)
axs[0].fill_between([1.5, 10.5], [5 - width, 5 - width], [5 + width, 5 + width], facecolor = 'dodgerblue', alpha = 0.3)


axs[0].fill_between([11.5, 18.5], [8 - width, 8 - width* 0.2], [8 + width, 8 + width * 0.2], facecolor = 'deeppink', alpha = 0.3)
axs[0].fill_between([11.5, 18.5], [2 - width, 2 - width* 0.2], [2 + width, 2 + width * 0.2], facecolor = 'darkorange', alpha = 0.3)
axs[0].fill_between([11.5, 18.5], [5 - width, 5 - width* 0.2], [5 + width, 5 + width * 0.2], facecolor = 'dodgerblue', alpha = 0.3)


#Etalon Area used

au1 = Rectangle((6 - 0.25, 4.33), 0.5, 1.33, facecolor = 'deeppink', alpha = 0.7  , edgecolor = 'white', linewidth =1)
au2 = Rectangle((6 - 0.25, 4.33), 0.5, 1.33, facecolor = 'darkorange', alpha = 0.7, edgecolor = 'white', linewidth =1)
au3 = Rectangle((6 - 0.25, 4.33), 0.5, 1.33, facecolor = 'dodgerblue', alpha = 0.7, edgecolor = 'white', linewidth =1)




axs[0].set_ylim(0, 10)
axs[0].set_xlim(0, 20)

axs[0].add_patch(source)
axs[0].add_patch(etalon)
axs[0].add_patch(lens)
axs[0].add_patch(sensor)
axs[0].add_patch(au1)
axs[0].add_patch(au2)
axs[0].add_patch(au3)



axs[0].annotate("Object\nplane", xy=(1, 0.15), ha = "center")
axs[0].annotate("Etalon\n(at pupil)", xy=(6, 0.15), ha = "center")
axs[0].annotate("Lens", xy=(11, 0.3), ha = "center")
axs[0].annotate("Sensor", xy=(19, 0.3), ha = "center")




 # -----------------------------# 
 # TELE


source = Rectangle((1 - 0.25, 1), 0.5, 8, facecolor = 'k')
etalon = Rectangle((16 - 0.25, 1), 0.5, 8, facecolor = 'k')
lens = Rectangle((11 - 0.25, 1), 0.5, 8, facecolor = 'k')
sensor = Rectangle((19 - 0.25, 1), 0.5, 8, facecolor = 'k')


lr1 = axs[1].annotate("", xy=(10.5, 8), xytext=(1.5, 2),arrowprops=dict(arrowstyle="->", color = 'k', alpha = 0.5))
lr2 = axs[1].annotate("", xy=(10.5, 5), xytext=(1.5, 5),arrowprops=dict(arrowstyle="->", color = 'k', alpha = 0.5))
lr3 = axs[1].annotate("", xy=(10.5, 2), xytext=(1.5, 8),arrowprops=dict(arrowstyle="->", color = 'k', alpha = 0.5))


flr1 = axs[1].annotate("", xy=(18.5, 8), xytext=(11.5, 8),arrowprops=dict(arrowstyle="->", color = 'k', alpha = 0.5))
flr2 = axs[1].annotate("", xy=(18.5, 5), xytext=(11.5, 5),arrowprops=dict(arrowstyle="->", color = 'k', alpha = 0.5))
flr3 = axs[1].annotate("", xy=(18.5, 2), xytext=(11.5, 2),arrowprops=dict(arrowstyle="->", color = 'k', alpha = 0.5))


width = 0.5
axs[1].fill_between([1.5, 10.5], [2 - width, 8 - width], [2 + width, 8 + width], facecolor = 'deeppink', alpha = 0.3)
axs[1].fill_between([1.5, 10.5], [8 - width, 2 - width], [8 + width, 2 + width], facecolor = 'darkorange', alpha = 0.3)
axs[1].fill_between([1.5, 10.5], [5 - width, 5 - width], [5 + width, 5 + width], facecolor = 'dodgerblue', alpha = 0.3)


axs[1].fill_between([11.5, 18.5], [8 - width, 8 - width* 0.2], [8 + width, 8 + width* 0.2], facecolor = 'deeppink', alpha = 0.3)
axs[1].fill_between([11.5, 18.5], [2 - width, 2 - width* 0.2], [2 + width, 2 + width* 0.2], facecolor = 'darkorange', alpha = 0.3)
axs[1].fill_between([11.5, 18.5], [5 - width, 5 - width* 0.2], [5 + width, 5 + width* 0.2], facecolor = 'dodgerblue', alpha = 0.3)


axs[1].plot([6, 6], [1, 9], ls = '--', c = 'k', lw = 2)

axs[1].set_ylim(0, 10)
axs[1].set_xlim(0, 20)

axs[1].add_patch(source)
axs[1].add_patch(etalon)
axs[1].add_patch(lens)
axs[1].add_patch(sensor)






au1 = Polygon([(15.75, 8.25), (16.25, 8.23), (16.25,7.77), (15.75, 7.742),],facecolor = 'deeppink',    alpha = 1,  edgecolor = 'white', linewidth =1.3)
au2 = Polygon([(15.75, 8.25-3), (16.25, 8.23-3), (16.25,7.77-3), (15.75, 7.742-3),],facecolor = 'dodgerblue',    alpha = 1,  edgecolor = 'white', linewidth =1.3)
au3 = Polygon([(15.75, 8.25-6), (16.25, 8.23-6), (16.25,7.77-6), (15.75, 7.742-6),],facecolor = 'darkorange',    alpha = 1,  edgecolor = 'white', linewidth =1.3)


#au1 = Rectangle((16 - 0.25, 8 - 0.5), 0.5, width * 2, facecolor = 'deeppink',    alpha = 1,  edgecolor = 'white', linewidth =1.3)
#au2 = Rectangle((16 - 0.25, 5 - 0.5), 0.5, width * 2, facecolor = 'dodgerblue',  alpha = 1,  edgecolor = 'white', linewidth =1.3)
#au3 = Rectangle((16 - 0.25, 2 - 0.5), 0.5, width * 2, facecolor = 'darkorange',  alpha = 1,  edgecolor = 'white', linewidth =1.3)

axs[1].add_patch(au1)
axs[1].add_patch(au2)
axs[1].add_patch(au3)


axs[1].annotate("Object\nplane", xy=(1, 0.15), ha = "center")
axs[1].annotate("Pupil", xy=(6, 0.3), ha = "center")
axs[1].annotate("Lens", xy=(11, 0.3), ha = "center")
axs[1].annotate("Etalon", xy=(16, 0.3), ha = "center")
axs[1].annotate("Sensor", xy=(19, 0.3), ha = "center")



axs[0].axis("off")
axs[1].axis("off")


axs[0].set_title("Collimated")
axs[1].set_title("Telecentric")


plt.tight_layout()

plt.savefig("Plots/plots/EtalonConfigurations.pdf", bbox_inches = 'tight')
