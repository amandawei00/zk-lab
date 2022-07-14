import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys

sys.path.append('../data')
plt.rcParams['lines.linewidth'] = 0.75

csfont = {'fontname':'Times New Roman'}

def equal_arr(a1, a2):
    a1.sort()
    a2.sort()

    for i in range(0, len(a1)):
        if a1[i] != a2[i]:
            # print(str(a1[i]) + ", " + str(a2[i]))
            return False

    return True

# plots ratios of two solutions of BK equation at specified rapidities
def ratios(df1, df2, rap):
    for i in range(len(rap)):
        sub_ = df1.loc[df1['y'] == rap[i]]
        sub2_ = df2.loc[df2['y'] == rap[i]]

        print(sub_)
        print(sub2_)
        vr = np.array(sub_[['vr']])
        vfr = np.array(sub_[['vfr']])

        vr2 = np.array(sub2_[['vr']])
        vfr2 = np.array(sub2_[['vfr']])

        rat = vfr/vfr2
        print(np.concatenate(rat))
        plt.plot(vr, rat, label='y='+str(rap[i]))

    plt.legend()
    plt.xlabel('r')
    plt.ylabel('MV(amanda)/MV(heikki)')
    plt.xscale('log')
    plt.show()
    return 0

# plots BK solution at specified rapdities
def plot_bk(df, rap):
    for i in range(len(rap)):
        sub_ = df.loc[(df['y'] == rap[i])]

        vr = np.array(sub_[['vr']])
        vfr = np.array(sub_[['vfr']])

        plt.plot(vr, vfr, label="y = " + str(rap[i]), color="black")

    plt.xlabel("r")
    plt.ylabel("N(r,Y)")
    plt.xscale('log')
    plt.legend()
    plt.show()

def get_data(df, rap):

    sub_ = df.loc[(df['y'] == rap)]
    vr = np.array(sub_[['vr']])
    vfr = np.array(sub_[['vfr']])

    return [vr, vfr]

def load_df(name):
    df = pd.read_csv(name, sep='\t', comment='#', header=None)
    # df = pd.read_csv(name, delim_whitespace=True)
    df.columns = ['y', 'vr', 'vfr']

    df['y'] = (df['y'].astype('float32')).round(decimals=2)
    df['vr'] = df['vr'].astype('float64')
    df['vfr'] = df['vfr'].astype('float64')

    return df

def make_animation(df):
    fig = plt.figure()
    ax = plt.axes(xlim=(3e-7, 100), ylim=(0.0, 1.2))
    ax.set_xscale('log')
    line, = ax.plot([], [], lw=2)

    ax.set_xlabel('r')
    ax.set_ylabel('N(r, Y)')
    
    def init():
        dat = get_data(df, 0.0)
        line.set_data(dat[0], dat[1])
        return line,

    def animate(i):
        dat = get_data(df, i)
        line.set_data(dat[0], dat[1])
        return line,

    anim = animation.FuncAnimation(fig, animate, frames=np.arange(0.0, 30.0, 0.3), init_func=init, blit=True, interval=0.5)
    plt.show()


if __name__ == '__main__':

    # df4 = load_df('../results/RK2/bk_MVg.csv')
    # df4 = load_df('../results/RK4/bk_MVg.csv')
    df4 = load_df('../results/RK4/bk_MVg_heikki.csv')
    df5 = load_df('../../../rcbk/mvg_test.csv')
    rap = [0.0, 5., 10., 15.]
    ratios(df4, df5, rap)

