drellYan__ = ["DrellYanLO"]
drellYanAll__ = ["DrellYanTauTauLO","DrellYanLO","DrellYanTauTauLOHT0to100","DrellYanLOHT0to100"]
drellYanTauTau__ = ["DrellYanTauTauLO"]
ttbar__ = ["TT_Powheg"]
other__ = ["SingleTop", "Rare", "Diboson", "DrellYanTauTauLO", "DrellYanTauTauLOHT0to100", "WJets"]
def getBackgrounds(*samples):
    bkgs = []
    for sample in samples:
        if sample == "DY":
            bkgs.extend(drellYan__)
        elif sample == "TT":
            bkgs.extend(ttbar__)
        elif sample == "DYAll":
            bkgs.extend(drellYanAll__)
        elif sample == "DYTauTau":
            bkgs.extend(drellYanTauTau__)
    return bkgs
