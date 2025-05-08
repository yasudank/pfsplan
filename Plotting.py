import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import datetime
from astropy.time import Time
import astropy.units as u
from collections import OrderedDict


# Define distinct colormaps for each wg
cmap_dict = {
    'GE': LinearSegmentedColormap.from_list('ge_cmap', ['red', 'orange']),
    'GA': LinearSegmentedColormap.from_list('ga_cmap', ['yellow', 'green']),
    'CO': LinearSegmentedColormap.from_list('co_cmap', ['blue', 'violet'])
}

def get_color_for_point(wg, i, length):
    if length != 1:
        frac = i / float(length-1)
    else:
        frac = 0.5
    return cmap_dict[wg](frac)

def plotSchedule(schedule, dates, obscond, targetList, observer, priority=-1):
    nrow = max(len(dates), 7)
    fig, axes = plt.subplots(nrow, 1, figsize=(10, 7), sharex=True, sharey=True)

    for k in range(len(dates)):
        xx = list()
        yy = list()
        zz = list()
        for slot in schedule.get_used_slots():
            if slot.date != dates[k]:
                continue
            t = slot.target
            tt = datetime.datetime.strptime((slot.mid+observer.utcoffset).iso[:-4], '%Y-%m-%d %H:%M:%S').time()
            _tt = tt.hour + tt.minute/60
            if _tt < 15:
                _tt += 24
            xx.append(_tt)
            yy.append(obscond.altaz(slot.index, t.name).alt.deg)
            zz.append(get_color_for_point(t.wg, targetList.wg_objects[t.wg].index(t.name), len(targetList.wg_objects[t.wg])))

        xm = list()
        ma = list()
        for slot in schedule:
            if slot.date != dates[k]:
                continue
            tt = datetime.datetime.strptime((slot.mid+observer.utcoffset).iso[:-4], '%Y-%m-%d %H:%M:%S').time()
            _tt = tt.hour + tt.minute/60
            if _tt < 15:
                _tt += 24
            xm.append(_tt)
            ma.append(obscond.moon_altaz(slot.index).alt.deg)
                
        axes[k].scatter(xx, yy, c=zz)
        axes[k].plot(xm, ma, c='gray', linestyle='--')

        time = Time(dates[k]+' 12:00:00') - observer.utcoffset
        for _horizon in [-18, -12, -6, 0]:
            sun_set = observer.sun_set_time(time, which='next', horizon=_horizon*u.degree)
            sun_rise = observer.sun_rise_time(time, which='next', horizon=_horizon*u.degree)
            _set = (sun_set+observer.utcoffset).datetime.hour + (sun_set+observer.utcoffset).datetime.minute/60
            _rise = (sun_rise+observer.utcoffset).datetime.hour + (sun_rise+observer.utcoffset).datetime.minute/60
            if _rise < 15:
                _rise += 24
            axes[k].axvspan(18, _set, color='gold', alpha=0.25)
            axes[k].axvspan(_rise, 30, color='gold', alpha=0.25)

        axes[k].set_yticks(range(0, 91, 30))
        axes[k].set_xlim(18, 30)
        axes[k].set_ylim(0, 90)
        axes[k].text(18.05, 97, dates[k], fontsize=10)
        axes[k].grid()
        if k == len(dates)-1:
            axes[k].set_xlabel('Time (HST)')

        axes[k].fill_between([18,30], 0, 30, color='red', alpha=0.3)

    for k in range(len(dates), nrow):
        axes[k].axis('off')

    fig.text(0.07, 0.5, 'Elevation (degree)', va='center', rotation='vertical')
    
    plt.subplots_adjust(hspace=0.4)

    if priority >= 0:
        plt.savefig(f'elevation_{priority}.png')
    else:
        plt.savefig(f'elevation.png')

def plotSchedule_rotang(schedule, dates, obscond, targetList, observer, priority=-1):
    nrow = max(len(dates), 7)
    fig, axes = plt.subplots(nrow, 1, figsize=(10, 7), sharex=True, sharey=True)

    for k in range(len(dates)):
        xx = list()
        ys = list()
        ye = list()
        zz = list()
        for slot in schedule.get_used_slots():
            if slot.date != dates[k]:
                continue
            t = slot.target
            tt = datetime.datetime.strptime((slot.mid+observer.utcoffset).iso[:-4], '%Y-%m-%d %H:%M:%S').time()
            _tt = tt.hour + tt.minute/60
            if _tt < 15:
                _tt += 24
            xx.append(_tt)
            ys.append(obscond.rotang_start(slot.index, t.name).degree)
            ye.append(obscond.rotang_end(slot.index, t.name).degree)
            zz.append(get_color_for_point(t.wg, targetList.wg_objects[t.wg].index(t.name), len(targetList.wg_objects[t.wg])))

        axes[k].scatter(xx, ys, c=zz)
        axes[k].scatter(xx, ye, c=zz)

        time = Time(dates[k]+' 12:00:00') - observer.utcoffset
        for _horizon in [-18, -12, -6, 0]:
            sun_set = observer.sun_set_time(time, which='next', horizon=_horizon*u.degree)
            sun_rise = observer.sun_rise_time(time, which='next', horizon=_horizon*u.degree)
            _set = (sun_set+observer.utcoffset).datetime.hour + (sun_set+observer.utcoffset).datetime.minute/60
            _rise = (sun_rise+observer.utcoffset).datetime.hour + (sun_rise+observer.utcoffset).datetime.minute/60
            if _rise < 15:
                _rise += 24
            axes[k].axvspan(18, _set, color='gold', alpha=0.25)
            axes[k].axvspan(_rise, 30, color='gold', alpha=0.25)

        axes[k].set_yticks(range(-180, 181, 90))
        axes[k].set_xlim(18, 30)
        axes[k].set_ylim(-180, 180)
        axes[k].text(18.05, 185, dates[k], fontsize=10)
        axes[k].grid()
        if k == len(dates)-1:
            axes[k].set_xlabel('Time (HST)')

        axes[k].fill_between([18,30], -180, -174, color='red', alpha=0.3)
        axes[k].fill_between([18,30],  174,  180, color='red', alpha=0.3)

    for k in range(len(dates), nrow):
        axes[k].axis('off')

    fig.text(0.07, 0.5, 'Rot Angle (degree)', va='center', rotation='vertical')
    
    plt.subplots_adjust(hspace=0.4)

    if priority >= 0:
        plt.savefig(f'rot_angle_{priority}.png')
    else:
        plt.savefig(f'rot_angle.png')

def plotSchedule_ha(schedule, dates, obscond, targetList, observer, priority=-1):
    nrow = max(len(dates), 7)
    fig, axes = plt.subplots(nrow, 1, figsize=(10, 7), sharex=True, sharey=True)

    for k in range(len(dates)):
        xx = list()
        yy = list()
        zz = list()
        for slot in schedule.get_used_slots():
            if slot.date != dates[k]:
                continue
            t = slot.target
            tt = datetime.datetime.strptime((slot.mid+observer.utcoffset).iso[:-4], '%Y-%m-%d %H:%M:%S').time()
            _tt = tt.hour + tt.minute/60
            if _tt < 15:
                _tt += 24
            xx.append(_tt)
            yy.append(obscond.ha(slot.index, t.name).value)
            zz.append(get_color_for_point(t.wg, targetList.wg_objects[t.wg].index(t.name), len(targetList.wg_objects[t.wg])))

        axes[k].scatter(xx, yy, c=zz)

        time = Time(dates[k]+' 12:00:00') - observer.utcoffset
        for _horizon in [-18, -12, -6, 0]:
            sun_set = observer.sun_set_time(time, which='next', horizon=_horizon*u.degree)
            sun_rise = observer.sun_rise_time(time, which='next', horizon=_horizon*u.degree)
            _set = (sun_set+observer.utcoffset).datetime.hour + (sun_set+observer.utcoffset).datetime.minute/60
            _rise = (sun_rise+observer.utcoffset).datetime.hour + (sun_rise+observer.utcoffset).datetime.minute/60
            if _rise < 15:
                _rise += 24
            axes[k].axvspan(18, _set, color='gold', alpha=0.25)
            axes[k].axvspan(_rise, 30, color='gold', alpha=0.25)

        axes[k].set_yticks(range(-4, 6, 2))
        axes[k].set_xlim(18, 30)
        axes[k].set_ylim(-5, 5)
        axes[k].text(18.05, 5.06, dates[k], fontsize=10)
        axes[k].grid()
        if k == len(dates)-1:
            axes[k].set_xlabel('Time (HST)')

        axes[k].fill_between([18,30], -0.3, 0.3, color='red', alpha=0.3)

    for k in range(len(dates), nrow):
        axes[k].axis('off')

    fig.text(0.07, 0.5, 'Hour Angle (hour)', va='center', rotation='vertical')
    
    plt.subplots_adjust(hspace=0.4)

    if priority >= 0:
        plt.savefig(f'hour_angle_{priority}.png')
    else:
        plt.savefig(f'hour_angle.png')

def plotObservedCounts(targetList):
    nexp_wg_finished = OrderedDict(sorted(targetList.nexp_wg_finished.items(), key=lambda x: x[0]))
    print(nexp_wg_finished)

    nexp_wg_frac = {w: nexp_wg_finished[w] / sum(nexp_wg_finished.values()) * 100 for w in nexp_wg_finished.keys()}
    print({k: f"{v:.1f}%" for k, v in nexp_wg_frac.items()})
    print(f"Total {sum(nexp_wg_finished.values())} exposures are finished")

    xx = list()
    yy = list()
    zz = list()
    for t in targetList.get_all_targets():
        if t.observed != 0 or t.wg != 'CO':
            xx.append(t.name)
            yy.append(t.observed)
            zz.append(get_color_for_point(t.wg, targetList.wg_objects[t.wg].index(t.name), len(targetList.wg_objects[t.wg])))

    # 1. xx のターゲット名をキー、元のインデックスを値とする辞書を作成
    #    これにより、ターゲット名からインデックスを O(1) (平均) で検索できます。
    name_to_original_index = {name: i for i, name in enumerate(xx)}

    # 2. targetList.wg_objects の順序に基づいて、並べ替え後のインデックスリストを作成
    sorted_indices = []
    for wg in targetList.wg_list:
        # targetList.wg_objects[wg] には、そのWGのターゲット名が望ましい順序で入っている想定
        for tname in targetList.wg_objects[wg]:
            # name_to_original_index を使って、元のインデックスを高速に取得
            original_index = name_to_original_index.get(tname)
            # もし tname が xx に存在すれば (つまり original_index が None でなければ)
            if original_index is not None:
                sorted_indices.append(original_index)
            # else:
            #   tname が xx に存在しない場合の処理 (必要であれば)
            #   例えば、エラーログを出力するなど
            #   print(f"Warning: Target '{tname}' from wg_objects not found in xx.")

    # 3. sorted_indices を使って、xx, yy, zz を並べ替えた新しいリストを作成
    #    リスト内包表記を使うと効率的です。
    xx_sorted = [xx[i] for i in sorted_indices]
    yy_sorted = [yy[i] for i in sorted_indices]
    zz_sorted = [zz[i] for i in sorted_indices]

    # 4. 元の変数名を置き換える (必要に応じて)
    xx = xx_sorted
    yy = yy_sorted
    zz = zz_sorted

    # Create a bar plot
    fig, ax = plt.subplots(figsize=(10, 7))
    bars = ax.bar(xx, yy, color=zz)
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
    # Set Y-axis limit with some padding
    max_height = max(yy) if yy else 0  # Handle empty data case
    ax.set_ylim(0, max_height * 1.1)  # Add 10% padding
    ax.set_xlabel('Target Name')
    ax.set_ylabel('Observed Counts')
    ax.set_title('Observed Counts for Each Target')
    plt.xticks(rotation=90)
    plt.tight_layout()

    wg_colors = {'GE': 'orangered', 'GA': 'yellowgreen', 'CO': 'blueviolet'}
    ax_inset = fig.add_axes([0.2, 0.51, 0.3, 0.35])
    bars_inset = ax_inset.bar(nexp_wg_finished.keys(), nexp_wg_finished.values(),
                             color=[wg_colors[w] for w in nexp_wg_finished.keys()])
    for i, bar in enumerate(bars_inset):
        height = bar.get_height()
        ax_inset.annotate(f'{height}',
                          xy=(bar.get_x() + bar.get_width() / 2, height),
                          xytext=(0, 3),  # 3 points vertical offset
                          textcoords="offset points",
                          ha='center', va='bottom')
        ax_inset.annotate(f'{nexp_wg_frac[list(nexp_wg_finished.keys())[i]]:.1f}%',
                            xy=(bar.get_x() + bar.get_width() / 2, 0),
                            xytext=(0, 3),  # 3 points vertical offset
                            textcoords="offset points",
                            ha='center', va='bottom')
    ax_inset.set_xlabel('Working Group')
    ax_inset.set_ylabel('Observed Counts')
    ax_inset.set_title('Observed Counts by WG')
    ax_inset.set_ylim(0, max(nexp_wg_finished.values()) * 1.1)

    plt.savefig('observed_counts.png')