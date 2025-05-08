from Params import Params
from MyObserver import MyObserver
from ObsSlot import ObsDate
from Targets import TargetManager, TargetList, Target
from ObservingConditions import ObservingConditions
from Optimize import OptimizeSchedule
from Plotting import plotSchedule, plotSchedule_rotang, plotSchedule_ha, plotObservedCounts
from Report import printSchedule as report_printSchedule, printSchedule_PDF # Renamed to avoid conflict
import logging
import pprint
import argparse # Added for command-line arguments

logger = logging.getLogger(__name__) # Logger can be defined globally


def optimization_1st(targetList, obsSlotList, ObservingConditions, params, subaru, obsdate):
    """
    1st stage of the optimization process.
    """
    for priority in targetList.priorities:
        logger.info(f"Optimization 1st stage - Priority: {priority}")
        o, obs_slots, targets, dummy = OptimizeSchedule(obsSlotList, targetList, ObservingConditions, params, subaru, obsdate.nexp_max, priority)

        #for slot in obs_slots:
        #    for t in targets + [dummy]:
        #        if o[slot.index, t.name].varValue > 0.5:
        #            print(f'{slot.date} {slot.index:3d} {t.name} {t.wg} {t.priority} {t.observed}')
        #            # logger.debug(f'Slot assignment (1st opt): {slot.date} {slot.index:3d} {t.name} {t.wg} {t.priority} {t.observed}')

        obsSlotList.updateSchedule(o, obs_slots, targets, targetList)
        #for t in targets:
        #    print(f'{t.name} {t.wg} {t.priority} {t.observed}')
        #    # logger.debug(f'Target observed (1st opt): {t.name} {t.wg} {t.priority} {t.observed}')

        plotSchedule(obsSlotList, obsdate.dates_local, ObservingConditions, targetList, subaru, priority)
        plotSchedule_rotang(obsSlotList, obsdate.dates_local, ObservingConditions, targetList, subaru, priority)
        plotSchedule_ha(obsSlotList, obsdate.dates_local, ObservingConditions, targetList, subaru, priority)

def optimization_2nd(obsSlotList, ObservingConditions, params, subaru, obsdate):
    """
    2nd stage of the optimization process.
    """
    # Create a new target list for the second stage
    targetList2 = TargetList()
    for slot in obsSlotList.get_used_slots():
        target = slot.target
        if not target.name in targetList2.names:
            target2 = Target(target.wg, target.name, target.coord, target.pa, target.observed, target.priority)
            targetList2.add_target(target2)

    logger.info(f'Total {targetList2.num_targets} targets available for 2nd stage')
    logger.info(f"WG objects for 2nd stage: {pprint.pformat(targetList2.wg_objects)}")

    # Load observation slots for the second stage
    obsdate2 = ObsDate(params.fname_obsdate,
                      params.fname_obsdate_finish,
                      observer=subaru,
                      params=params)
    obsSlotList2 = obsdate2.obsSlotList

    o, obs_slots, targets, dummy = OptimizeSchedule(obsSlotList2, targetList2, ObservingConditions, params, subaru, obsdate.nexp_max)
    obsSlotList2.updateSchedule(o, obs_slots, targets, targetList2)

    logger.info("Targets after 2nd optimization stage:")
    for t in targets:
        logger.info(f"  Name: {t.name}, WG: {t.wg}, Priority: {t.priority}, Observed: {t.observed}")

    plotSchedule(obsSlotList2, obsdate.dates_local, ObservingConditions, targetList2, subaru)
    plotSchedule_rotang(obsSlotList2, obsdate.dates_local, ObservingConditions, targetList2, subaru)
    plotSchedule_ha(obsSlotList2, obsdate.dates_local, ObservingConditions, targetList2, subaru)

    return obsSlotList2, targetList2

def optimization_3rd(obsSlotList, ObservingConditions, params, subaru, obsdate):
    """
    2nd stage of the optimization process.
    """
    # Create a new target list for the second stage
    targetList = TargetList()
    for slot in obsSlotList.get_used_slots():
        target = slot.target
        if not target.name in targetList.names:
            target = Target(target.wg, target.name, target.coord, target.pa, target.observed, target.priority)
            targetList.add_target(target)

    logger.info(f'Total {targetList.num_targets} targets available for 3rd stage')
    logger.info(f"WG objects for 3rd stage: {pprint.pformat(targetList.wg_objects)}")

    obsSlotList.reset()

    o, obs_slots, targets, dummy = OptimizeSchedule(obsSlotList, targetList, ObservingConditions, params, subaru, obsdate.nexp_max)
    obsSlotList.updateSchedule(o, obs_slots, targets, targetList)

    logger.info("Targets after 3rd optimization stage:")
    for t in targets:
        logger.info(f"  Name: {t.name}, WG: {t.wg}, Priority: {t.priority}, Observed: {t.observed}")

    plotSchedule(obsSlotList, obsdate.dates_local, ObservingConditions, targetList, subaru)
    plotSchedule_rotang(obsSlotList, obsdate.dates_local, ObservingConditions, targetList, subaru)
    plotSchedule_ha(obsSlotList, obsdate.dates_local, ObservingConditions, targetList, subaru)

    return obsSlotList, targetList

def reorderGAtargets(obsSlotList):
    GA_fields = list(set([slot.target.name [:slot.target.name.rindex('_')] for slot in obsSlotList.get_used_slots() if slot.target.wg == 'GA']))
    #logger.debug(f"GA fields for reordering: {GA_fields}")
    for ga_field in GA_fields:
        #logger.debug(f"Reordering GA field: {ga_field}")
        ga_slots = obsSlotList.get_slots_by_field(ga_field)
        ga_targets = sorted(list(set([slot.target for slot in ga_slots])), key=lambda x: (x.priority, x.name))
        #logger.debug(f"GA field {ga_field}, sorted targets: {[t.name for t in ga_targets]}")
        #logger.debug(f"GA field {ga_field}, number of slots: {len(ga_slots)}")
        i0 = 0
        for t in ga_targets:
            nobs = 0
            for i in range(t.nexp):
                if i0 + i >= len(ga_slots):
                    break
                ga_slots[i0 + i].target = t
                nobs += 1
            t.observed = nobs
            i0 += nobs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="SSP Plan Optimizer")
    parser.add_argument(
        '--log-level',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        help='Set the logging level for the application (default: INFO)'
    )
    args = parser.parse_args()

    # Configure logging using the command-line argument
    log_level_numeric = getattr(logging, args.log_level.upper(), None)
    if not isinstance(log_level_numeric, int):
        # This should not happen with choices, but as a safeguard
        raise ValueError(f'Invalid log level: {args.log_level}')

    logging.basicConfig(
        level=log_level_numeric,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    # Load parameters
    params = Params('parameters_2025May.yaml')

    logger.info(f"Airmass limits: {params.airmass['limit']}")

    # Initialize the observer
    subaru = MyObserver.at_site('Subaru', timezone='US/Hawaii')
    logger.info(f"UTC offset: {subaru.utcoffset}")

    # Load observation slots
    obsdate = ObsDate(params.fname_obsdate,
                      params.fname_obsdate_finish,
                      observer=subaru,
                      params=params)
    obsSlotList = obsdate.obsSlotList
    num_slots = obsSlotList.num_slots
    logger.info(f"Total {num_slots} observation slots available")
    
    for date in obsdate.dates_local:
        logger.info(f"Slots for date {date}: {[slot.index for slot in obsdate.obsSlotList.get_slots_by_date(date)]}")

    logger.info(f"Nexp max per WG: {obsdate.nexp_max}")
    logger.info(f"Nexp min per WG: {obsdate.nexp_min}")
    
    # Load target list
    target_manager = TargetManager(params.fname_targets, params.fname_targets_finish)
    targetList = target_manager.targetList
    num_targets = targetList.num_targets
    logger.info(f'Total {num_targets} targets available')
    logger.info(f"Priorities: {targetList.priorities}")
    logger.info(f"Working groups: {targetList.wg_list}")
    logger.info(f"WG objects: {pprint.pformat(targetList.wg_objects)}")
    #logger.debug(pprint.pformat(targetList.wg_objects)) # if pprint.pprint was used for debugging


    observingConditions = ObservingConditions(obsSlotList, targetList, subaru, params)

    # 1st stage of the optimization
    optimization_1st(targetList, obsSlotList, observingConditions, params, subaru, obsdate)

    reorderGAtargets(obsSlotList)

    # 2nd stage of the optimization
    obsSlotList2, targetList2 = optimization_2nd(obsSlotList, observingConditions, params, subaru, obsdate)
    
    for t in targetList.get_all_targets():
        if not t.name in targetList2.names:
            targetList2.add_target(t)

    reorderGAtargets(obsSlotList2)

    plotObservedCounts(targetList2)

    report_printSchedule(obsSlotList2, obsdate.dates_local, subaru, observingConditions, params)

    obsSlotList2.updateTimeBySlew(observingConditions, params)

    observingConditions2 = ObservingConditions(obsSlotList2, targetList2, subaru, params) # Recalculate OC with updated times

    # 3rd stage of the optimization
    #obsSlotList3, targetList3 = optimization_3rd(obsSlotList2, observingConditions2, params, subaru, obsdate)

    #for t in targetList.get_all_targets():
    #    if not t.name in targetList3.names:
    #        targetList3.add_target(t)

    report_printSchedule(obsSlotList2, obsdate.dates_local, subaru, observingConditions2, params)

    import pickle
    with open('obsSlotList2.pickle', 'wb') as f:
        pickle.dump(obsSlotList2, f)
    with open('observingConditions2.pickle', 'wb') as f:
        pickle.dump(observingConditions2, f)

    printSchedule_PDF(obsSlotList2, subaru, observingConditions2, params)