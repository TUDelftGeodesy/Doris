# This files contains a number of different network definitions, based temporal baselines, perpendicular baselines and
# a reference date (not always needed)

from datetime import datetime


def cascade_network(dates):
    # Sort based on datetimes and assign different

    dates = [datetime.strftime(d, '%Y-%m-%d') for d in dates]
    dates = sorted(dates)

    ifgs_pairs = []
    ifgs_keys = []
    for t1, t2 in zip(dates[:-1], dates[1:]):
        ifgs_pairs.append((t1, t2))
        ifgs_keys.append((datetime.strftime(t1,'%Y%m%d') + '_' + datetime.strftime(t2,'%Y%m%d')))

    return ifgs_pairs, ifgs_keys


def threshold_network(dates, perp_b, base_t, date_t):
    # Sort based on datetimes and assign different (perpendicular baseline in m and temporal baseline in days...)

    dates = [datetime.strptime(d, '%Y-%m-%d') for d in dates]
    ids = [i[0] for i in sorted(enumerate(dates), key=lambda x:x[1])]

    dates = [dates[idn] for idn in ids]
    perp_b = [perp_b[idn] for idn in ids]

    ifgs_pairs = []
    ifgs_keys = []
    for i, t1, d1 in zip(range(len(dates)), dates, perp_b):
        for t2, d2 in zip(dates[i+1:], perp_b[i+1:]):
            if abs((t2-t1).days) <= date_t and abs(d2-d1) <= base_t:
                ifgs_pairs.append((t1, t2))
                ifgs_keys.append(datetime.strftime(t1,'%Y%m%d') + '_' + datetime.strftime(t2,'%Y%m%d'))

    return ifgs_pairs, ifgs_keys


def exp_coh_network(dates, perp_b, coh_t, par_t, par_p):
    # Now the network is calculated based on coherence,

    print('in development')


def seasonal_coh_network(dates, perp_b, coh_t, ref_date):

    print('in development')


def exp_coherence():

    print('in development')


def seasonal_coherence():

    print('in development')
