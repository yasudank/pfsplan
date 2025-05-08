from pulp import LpVariable, LpProblem, LpMaximize, lpSum, LpStatus, PULP_CBC_CMD
from Targets import Target
import pprint
import logging

# Configure logging
logger = logging.getLogger(__name__)

def split_into_continuous_sequences(numbers):
    """
    Split a list of integers into sublists of continuous sequences.
    
    Args:
        numbers: List of integers
        
    Returns:
        List of lists, where each sublist contains a continuous sequence
    """
    if not numbers:
        return []
        
    # Sort the list first
    sorted_numbers = sorted(numbers)
    
    result = []
    current_sequence = [sorted_numbers[0]]
    
    # Iterate through the sorted list starting from the second number
    for i in range(1, len(sorted_numbers)):
        # If current number is consecutive to the previous one
        if sorted_numbers[i] == sorted_numbers[i-1] + 1:
            current_sequence.append(sorted_numbers[i])
        else:
            # End of a continuous sequence, start a new one
            result.append(current_sequence)
            current_sequence = [sorted_numbers[i]]
    
    # Don't forget to add the last sequence
    if current_sequence:
        result.append(current_sequence)
        
    return result

# Define the optimazation problem
def OptimizeSchedule(obsSlotList, targetList, oc, params, observer, nexp_max, priority=-1):

    dummy = Target('dummy', 'dummy', None, None, obsSlotList.num_slots, 10, 0)

    if priority == -1:
        targets = targetList.get_all_targets()
    else:        
        targets = targetList.get_observing_targets_by_priority(priority)
    logger.info(f'Total {len(targets)} targets available for priority targets (priority <= {priority})')
    for wg in targetList.wg_list:
        logger.info(f"  WG {wg}: {[t.name for t in targets if t.wg == wg]}")
    targets_with_dummy = targets + [dummy]
    wg_list_with_dummy = targetList.wg_list + [dummy.wg]

    nslot_required = sum([(t.nexp - t.observed) for t in targets])
    logger.info(f'Total {nslot_required} slots required for priority targets (priority <= {priority})')

    obs_slots = obsSlotList.get_available_slots(nslot_required)
    logger.info(f'Total {len(obs_slots)} slots available for priority targets (priority <= {priority})')
    
    dates_utc = list(set([slot.start.strftime('%Y-%m-%d') for slot in obs_slots]))
    dates_utc.sort()
    logger.info(f"Dates (UTC) for optimization: {dates_utc}")
    slots_by_date = {date_utc: [slot.index for slot in obs_slots if slot.start.strftime('%Y-%m-%d') == date_utc] for date_utc in dates_utc}
    #logger.debug(f"Slots by date: {pprint.pformat(slots_by_date)}")

    # Define the problem
    prob = LpProblem("ObservingPlan", LpMaximize)

    # Define the variables : o[slot, target] = 1 if the target is observed in the slot
    o = LpVariable.dicts('o', [(slot.index, t.name) for slot in obs_slots for t in targets_with_dummy], cat='Binary')

    # Define the Variable : y[date, target] = 1 if the obsevation of the target is completed on the date
    y = LpVariable.dicts('y', [(date, t.name) for date in dates_utc for t in targets_with_dummy], cat='Binary')

    # Objective function: maximize the sum of the effective exposure time of the targets that are observed
    prob += lpSum([o[(slot.index, t.name)] * oc.teff(slot.index, t.name) for slot in obs_slots for t in targets_with_dummy]) \
            + params.weight_comp * lpSum([y[(date, t.name)] for date in dates_utc for t in targets_with_dummy]) \
            - params.weight_pri * lpSum([t.priority * o[(slot.index, t.name)] for slot in obs_slots for t in targets]) \
            #- weight_slew * lpSum([x[(obs_slots[j].index, t1.name, t2.name)] * slewTime[(obs_slots[j].index, t1.name, t2.name)].to(u.minute).value for j in range(len(obs_slots)-1) if obs_slots[j]['date'] == obs_slots[j+1]['date'] for t1 in targets_with_dummy for t2 in targets_with_dummy])  # - slew time
            #- weight_slew * lpSum([o[(obs_slots[j].index, t1.name)] * o[(obs_slots[j+1].index, t2.name)] * slewTime[(obs_slots[j].index, t1.name, t2.name)].to(u.minute).value for j in range(len(obs_slots)-1) if obs_slots[j]['date'] == obs_slots[j+1]['date'] for t1 in targets_with_dummy for t2 in targets_with_dummy])  # - slew time

    # Constraints: each target is observed at most nexp times
    for t in targets_with_dummy:
        prob += lpSum([o[(slot.index, t.name)] for slot in obs_slots]) + t.observed <= t.nexp

    # Constraints: y = 1 if the target is observed nexp times on the date
    for date in dates_utc:
        for t in targets_with_dummy:
            prob += lpSum([o[(islot, t.name)] for islot in slots_by_date[date]]) >= t.nexp * y[(date, t.name)]
                      
    # Constraints: each timeslot is used for at most one target
    for slot in obs_slots:
        prob += lpSum([o[(slot.index, t.name)] for t in targets_with_dummy]) == 1

    # Constraints: Limitation for the number of exposures for each wg (Maximum)
    for w in targetList.wg_list:
        prob += lpSum([o[(slot.index, t.name)] for t in targets if t.wg == w for slot in obs_slots]) \
            <= nexp_max[w] - targetList.nexp_wg_finished[w]
    
    # Variable indicating whether each WG targets are observed (True) or not (False) at each timeslot
    wg_obs = LpVariable.dicts('wg_obs', [(date, i, w) for date in dates_utc for i in slots_by_date[date] for w in wg_list_with_dummy], cat='Bianry')
    for date in dates_utc:
        for i in slots_by_date[date]:
            for w in wg_list_with_dummy:
                prob += wg_obs[(date, i, w)] == lpSum([o[(i, t.name)] for t in targets_with_dummy if t.wg == w])
    
    
    # Constraints: GA targets will be observed after GE and CO targets
    if params.GA_last:
        for date in dates_utc:
            for i1 in slots_by_date[date]:
                for i2 in slots_by_date[date]:
                    if i1 > i2:
                        prob += wg_obs[(date, i1, "GA")] >= wg_obs[(date, i2, "GA")]
    
    
    # Variable indicating when the observation starts and ends
    wg_start = LpVariable.dicts('wg_start', [(date, i, w) for date in dates_utc for i in slots_by_date[date] for w in wg_list_with_dummy], \
                                lowBound=0, cat='Bianry')
    wg_end   = LpVariable.dicts('wg_end',   [(date, i, w) for date in dates_utc for i in slots_by_date[date] for w in wg_list_with_dummy], \
                                lowBound=0, cat='Bianry')
    
    # Constraints: Define wg_start and wg_end
    # Referene: https://techblog.zozo.com/entry/mip-wfm-scheduling
    for date in dates_utc:
        sub_slots = split_into_continuous_sequences(slots_by_date[date])
        for sub_slot in sub_slots:
            i_s = sub_slot[0]
            i_e = sub_slot[-1]
            for w in wg_list_with_dummy:
                prob += wg_obs[(date, i_s, w)] == wg_start[(date, i_s, w)]
                prob += wg_obs[(date, i_e, w)] == wg_end[(date, i_e, w)]
            if len(sub_slot) == 1:
                continue
            for i in sub_slot[1:]:
                for w in wg_list_with_dummy:
                    prob += wg_obs[(date, i, w)] - wg_obs[(date, i-1, w)] <= wg_start[(date, i, w)]
                    prob += wg_start[(date, i, w)] <= wg_obs[(date, i, w)]
                    prob += wg_start[(date, i, w)] <= 1 - wg_obs[(date, i-1, w)]
            for i in sub_slot[:-1]:
                for w in wg_list_with_dummy:
                    prob += wg_obs[(date, i, w)] - wg_obs[(date, i+1, w)] <= wg_end[(date, i, w)]
                    prob += wg_end[(date, i, w)] <= wg_obs[(date, i, w)]
                    prob += wg_end[(date, i, w)] <= 1 - wg_obs[(date, i+1, w)]
    
    
    # Constraints: The observation of each WG should continue at least n_continuous time slots
    """
    for date in dates_utc:
        sub_slots = split_into_continuous_sequences(slots_by_date[date])
        for sub_slot in sub_slots:
            if len(sub_slot) < params.n_continuous:
                logger.debug(f"Skipping continuous constraint for date {date}, sub_slot {sub_slot} (len {len(sub_slot)}) due to n_continuous {params.n_continuous}")
                continue
            for i in sub_slot[params.n_continuous-1:]:
                for w in wg_list_with_dummy:
                    for j in range(params.n_continuous):
                        prob += wg_end[(date, i, w)] <= wg_obs[(date, i-j, w)]
            for i in sub_slot[:-params.n_continuous+1]:
                for w in wg_list_with_dummy:
                    for j in range(params.n_continuous):
                        prob += wg_start[(date, i, w)] <= wg_obs[(date, i+j, w)]
    """

    # Constraints: The observation of each WG should continue at least n_continuous time slots
    for date in dates_utc:
        sub_slots = split_into_continuous_sequences(slots_by_date[date])
        for sub_slot in sub_slots:
            for w in wg_list_with_dummy:
                if len(sub_slot) < params.n_continuous[w]:
                    logger.debug(f"Skipping continuous constraint for WG {w} on date {date}, sub_slot {sub_slot} (len {len(sub_slot)}) due to n_continuous[{w}] {params.n_continuous[w]}")
                    continue
                for i in sub_slot[params.n_continuous[w]-1:]:
                    for j in range(params.n_continuous[w]):
                        prob += wg_end[(date, i, w)] <= wg_obs[(date, i-j, w)]
                for i in sub_slot[:-params.n_continuous[w]+1]:
                    for j in range(params.n_continuous[w]):
                        prob += wg_start[(date, i, w)] <= wg_obs[(date, i+j, w)]
    
    ######################################################################################################################
    if priority == -1:
        # Variable indicating whether each targets are observed (True) or not (False) at each timeslot
        tg_obs = LpVariable.dicts('tg_obs', [(date, i, t.name) for date in dates_utc for i in slots_by_date[date] for t in targets_with_dummy], cat='Bianry')
        for date in dates_utc:
            for i in slots_by_date[date]:
                for tg in targets_with_dummy:
                    prob += tg_obs[(date, i, tg.name)] == lpSum([o[(i, t.name)] for t in targets_with_dummy if t.name == tg.name])

        # Variable indicating when the observation starts and ends
        tg_start = LpVariable.dicts('tg_start', [(date, i, t.name) for date in dates_utc for i in slots_by_date[date] for t in targets_with_dummy], \
                                    lowBound=0, cat='Bianry')
        tg_end   = LpVariable.dicts('tg_end',   [(date, i, t.name) for date in dates_utc for i in slots_by_date[date] for t in targets_with_dummy], \
                                    lowBound=0, cat='Bianry')

        # Constraints: Define wg_start and wg_end
        # Referene: https://techblog.zozo.com/entry/mip-wfm-scheduling
        for date in dates_utc:
            sub_slots = split_into_continuous_sequences(slots_by_date[date])
            for sub_slot in sub_slots:
                i_s = sub_slot[0]
                i_e = sub_slot[-1]
                for t in targets_with_dummy:
                    prob += tg_obs[(date, i_s, t.name)] == tg_start[(date, i_s, t.name)]
                    prob += tg_obs[(date, i_e, t.name)] == tg_end[(date, i_e, t.name)]
                if len(sub_slot) == 1:
                    continue
                for i in sub_slot[1:]:
                    for t in targets_with_dummy:
                        prob += tg_obs[(date, i, t.name)] - tg_obs[(date, i-1, t.name)] <= tg_start[(date, i, t.name)]
                        prob += tg_start[(date, i, t.name)] <= tg_obs[(date, i, t.name)]
                        prob += tg_start[(date, i, t.name)] <= 1 - tg_obs[(date, i-1, t.name)]
                for i in sub_slot[:-1]:
                    for t in targets_with_dummy:
                        prob += tg_obs[(date, i, t.name)] - tg_obs[(date, i+1, t.name)] <= tg_end[(date, i, t.name)]
                        prob += tg_end[(date, i, t.name)] <= tg_obs[(date, i, t.name)]
                        prob += tg_end[(date, i, t.name)] <= 1 - tg_obs[(date, i+1, t.name)]

        # Constraints: The observation of each WG should continue at least n_continuous time slots
        n_cont_target = 2
        for date in dates_utc:
            sub_slots = split_into_continuous_sequences(slots_by_date[date])
            for sub_slot in sub_slots:
                if len(sub_slot) < n_cont_target:
                    logger.debug(f"Skipping target continuous constraint on date {date}, sub_slot {sub_slot} (len {len(sub_slot)}) due to n_cont_target {n_cont_target}")
                    continue
                for i in sub_slot[n_cont_target-1:]:
                    for t in targets_with_dummy:
                        for j in range(n_cont_target):
                            prob += tg_end[(date, i, t.name)] <= tg_obs[(date, i-j, t.name)]
                for i in sub_slot[:-n_cont_target+1]:
                    for t in targets_with_dummy:
                        for j in range(n_cont_target):
                            prob += tg_start[(date, i, t.name)] <= tg_obs[(date, i+j, t.name)]
    #
    ######################################################################################################################

    rotang_min = {'CO': -164, 'GE': -164, 'GA': -164}
    rotang_max = {'CO':  164, 'GE':  164, 'GA':  164}
    for slot in obs_slots:
        for t in targets:
            # Constraints: the moon is at least 60 degrees away from each target
            prob += oc.moon_sep(slot.index, t.name).degree * o[(slot.index, t.name)] \
                >= params.moonsep['limit'].value * o[(slot.index, t.name)]

            # Constraints: the Mars is at least 10 degrees away from each target
            for planet in ["mars", "jupiter", "saturn"]:
                prob += oc.planet_sep(planet, slot.index, t.name).degree * o[(slot.index, t.name)] \
                    >= params.planetssep['limit'].value * o[(slot.index, t.name)]
            
            # Constraints: the airmass is less than the limit
            prob += oc.airmass(slot.index, t.name).value * o[(slot.index, t.name)] \
                <= params.airmass['limit'][t.wg] * o[(slot.index, t.name)]

            # Constraints: we need to avoid the meridian
            prob += abs(oc.ha(slot.index, t.name).value) * o[(slot.index, t.name)] \
                >= params.meridian['warn'].value * o[(slot.index, t.name)] * (t.coord.dec > observer.location.lat)

            # Constraints: rotator limit
            prob += oc.rotang_start(slot.index, t.name).degree * o[(slot.index, t.name)] \
                >= rotang_min[t.wg] * o[(slot.index, t.name)]
            prob += oc.rotang_start(slot.index, t.name).degree * o[(slot.index, t.name)] \
                <= rotang_max[t.wg] * o[(slot.index, t.name)]
            prob += oc.rotang_end(slot.index, t.name).degree * o[(slot.index, t.name)] \
                >= rotang_min[t.wg] * o[(slot.index, t.name)]
            prob += oc.rotang_end(slot.index, t.name).degree * o[(slot.index, t.name)] \
                <= rotang_max[t.wg] * o[(slot.index, t.name)]
            
    # Solve the problem
    prob.solve(PULP_CBC_CMD(msg=0, threads=8))

    logger.info(f"Optimization status: {LpStatus[prob.status]}")

    return o, obs_slots, targets, dummy
