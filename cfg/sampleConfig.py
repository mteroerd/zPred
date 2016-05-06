drellYan__ = ["DrellYanLO", "DrellYanLOHT0to100"]
ttbar__ = ["TT_Powheg"]

def getBackgrounds(*samples):
    bkgs = []
    for sample in samples:
        if sample == "DY":
            bkgs.extend(drellYan__)
        elif sample == "TT":
            bkgs.extend(ttbar__)
    return bkgs
