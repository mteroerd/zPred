drellYan__ = ["DrellYanLO", "DrellYanLOHT0to100"]
drellYanTauTau__ = ["DrellYanTauTauLO", "DrellYanTauTauLOHT0to100"]
ttbar__ = ["TT_Powheg"]
other__ = ["SingleTop", "Rare", "Diboson", "DrellYanTauTauLO", "DrellYanTauTauLOHT0to100", "WJets"]
def getBackgrounds(*samples):
    bkgs = []
    for sample in samples:
        if sample == "DY":
            bkgs.extend(drellYan__)
        elif sample == "TT":
            bkgs.extend(ttbar__)
        elif sample == "DYTauTau":
            bkgs.extend(drellYanTauTau__)
    return bkgs
