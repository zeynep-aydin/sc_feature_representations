# fixed_splits.R
#
# Pre-defined leave-cohort-out train/test splits for datasets with few,
# heavily imbalanced donors where random splitting produces unstable fractions.
#
# Applies to: zilionis_lung, scea
# Does NOT apply to: he_organs, zhao_immune (no donor/batch metadata available)
#
# 10 splits per target fraction (50%, 70%, 80%), selected to:
#   - Land within ±2% of the target (±3% for Zilionis 80%, the limit of 7 donors)
#   - Maximise donor-role diversity across the 10 splits

FIXED_SPLITS <- list(

  # ── Zilionis Lung ────────────────────────────────────────────────────────
  # 7 donors: p1=28.1%, p3=20.7%, p7=12.5%, p6=11.0%, p4=9.9%, p2=9.4%, p5=8.4%
  zilionis_lung = list(
    "50%" = list(
      S01 = list(train=c("p3","p6","p4","p2"),         test=c("p1","p7","p5")),          # 51.0%
      S02 = list(train=c("p1","p7","p5"),               test=c("p3","p6","p4","p2")),     # 49.0%
      S03 = list(train=c("p1","p3"),                    test=c("p7","p6","p4","p2","p5")),# 48.9%
      S04 = list(train=c("p7","p6","p4","p2","p5"),     test=c("p1","p3")),               # 51.1%
      S05 = list(train=c("p1","p7","p6"),               test=c("p3","p4","p2","p5")),     # 51.6%
      S06 = list(train=c("p3","p4","p2","p5"),          test=c("p1","p7","p6")),          # 48.4%
      S07 = list(train=c("p1","p7","p4"),               test=c("p3","p6","p2","p5")),     # 50.5%
      S08 = list(train=c("p3","p6","p2","p5"),          test=c("p1","p7","p4")),          # 49.5%
      S09 = list(train=c("p1","p6","p4"),               test=c("p3","p7","p2","p5")),     # 49.0%
      S10 = list(train=c("p3","p7","p2","p5"),          test=c("p1","p6","p4"))           # 51.0%
    ),
    "70%" = list(
      S01 = list(train=c("p1","p3","p7","p4"),          test=c("p6","p2","p5")),          # 71.2%
      S02 = list(train=c("p1","p7","p6","p2","p5"),     test=c("p3","p4")),               # 69.4%
      S03 = list(train=c("p1","p3","p6","p4"),          test=c("p7","p2","p5")),          # 69.7%
      S04 = list(train=c("p3","p7","p6","p4","p2","p5"),test=c("p1")),                   # 71.9%
      S05 = list(train=c("p1","p3","p7","p2"),          test=c("p6","p4","p5")),          # 70.8%
      S06 = list(train=c("p1","p3","p6","p5"),          test=c("p7","p4","p2")),          # 68.2%
      S07 = list(train=c("p1","p7","p4","p2","p5"),     test=c("p3","p6")),               # 68.3%
      S08 = list(train=c("p1","p3","p6","p2"),          test=c("p7","p4","p5")),          # 69.3%
      S09 = list(train=c("p1","p7","p6","p4","p5"),     test=c("p3","p2")),               # 69.9%
      S10 = list(train=c("p1","p3","p4","p2"),          test=c("p7","p6","p5"))           # 68.2%
    ),
    "80%" = list(
      S01 = list(train=c("p1","p3","p7","p6","p4"),     test=c("p2","p5")),               # 82.2%
      S02 = list(train=c("p1","p3","p7","p2","p5"),     test=c("p6","p4")),               # 79.2%
      S03 = list(train=c("p1","p3","p6","p4","p2"),     test=c("p7","p5")),               # 79.1%
      S04 = list(train=c("p1","p7","p6","p4","p2","p5"),test=c("p3")),                   # 79.3%
      S05 = list(train=c("p1","p3","p7","p6","p5"),     test=c("p4","p2")),               # 80.7%
      S06 = list(train=c("p1","p3","p7","p4","p2"),     test=c("p6","p5")),               # 80.7%
      S07 = list(train=c("p1","p3","p6","p4","p5"),     test=c("p7","p2")),               # 78.0%
      S08 = list(train=c("p1","p3","p6","p2","p5"),     test=c("p7","p4")),               # 77.6%
      S09 = list(train=c("p1","p3","p7","p4","p5"),     test=c("p6","p2")),               # 79.6%
      S10 = list(train=c("p1","p3","p7","p6","p2"),     test=c("p4","p5"))                # 81.8%
    )
  ),

  # ── SCEA ─────────────────────────────────────────────────────────────────
  # 25 cohorts; E-ANND-1 dominates at 40.4%.
  # For 70% and 80%, E-ANND-1 must always be in train (mathematical constraint).
  scea = list(
    "50%" = list(
      S01 = list(
        train = c("E-HCAD-15","E-ANND-2","E-HCAD-1","E-HCAD-35","E-CURD-119","E-CURD-88",
                  "E-MTAB-10553","E-HCAD-31","E-GEOD-114530","E-HCAD-8","E-GEOD-125970","E-GEOD-130148"),
        test  = c("E-ANND-1","E-HCAD-10","E-MTAB-8410","E-CURD-122","E-HCAD-11","E-HCAD-9",
                  "E-MTAB-9543","E-GEOD-84465","E-ANND-3","E-GEOD-124472","E-MTAB-5061",
                  "E-MTAB-9067","E-GEOD-83139")),                                          # 50.0%
      S02 = list(
        train = c("E-ANND-1","E-CURD-119","E-HCAD-31","E-MTAB-8410","E-HCAD-11","E-HCAD-9",
                  "E-MTAB-9543","E-GEOD-84465","E-ANND-3","E-GEOD-124472","E-MTAB-5061",
                  "E-MTAB-9067","E-GEOD-83139"),
        test  = c("E-HCAD-15","E-ANND-2","E-HCAD-1","E-HCAD-35","E-HCAD-10","E-CURD-88",
                  "E-MTAB-10553","E-GEOD-114530","E-CURD-122","E-HCAD-8","E-GEOD-125970",
                  "E-GEOD-130148")),                                                        # 49.8%
      S03 = list(
        train = c("E-HCAD-15","E-ANND-2","E-HCAD-1","E-HCAD-35","E-HCAD-10","E-CURD-88",
                  "E-MTAB-10553","E-GEOD-114530","E-CURD-122","E-HCAD-8","E-GEOD-130148",
                  "E-ANND-3","E-GEOD-124472","E-MTAB-5061"),
        test  = c("E-ANND-1","E-CURD-119","E-HCAD-31","E-MTAB-8410","E-HCAD-11","E-HCAD-9",
                  "E-MTAB-9543","E-GEOD-125970","E-GEOD-84465","E-MTAB-9067","E-GEOD-83139")), # 50.5%
      S04 = list(
        train = c("E-ANND-1","E-HCAD-10","E-MTAB-8410","E-HCAD-11","E-HCAD-9","E-MTAB-9543",
                  "E-GEOD-125970","E-GEOD-84465","E-MTAB-5061","E-GEOD-83139"),
        test  = c("E-HCAD-15","E-ANND-2","E-HCAD-1","E-HCAD-35","E-CURD-119","E-CURD-88",
                  "E-MTAB-10553","E-HCAD-31","E-GEOD-114530","E-CURD-122","E-HCAD-8",
                  "E-GEOD-130148","E-ANND-3","E-GEOD-124472","E-MTAB-9067")),               # 48.5%
      S05 = list(
        train = c("E-HCAD-15","E-ANND-2","E-HCAD-1","E-HCAD-35","E-CURD-119","E-MTAB-10553",
                  "E-HCAD-31","E-GEOD-114530","E-MTAB-8410","E-CURD-122","E-HCAD-11",
                  "E-MTAB-9543","E-GEOD-84465","E-GEOD-130148","E-ANND-3","E-GEOD-124472",
                  "E-MTAB-9067","E-GEOD-83139"),
        test  = c("E-ANND-1","E-HCAD-10","E-CURD-88","E-HCAD-8","E-HCAD-9","E-GEOD-125970",
                  "E-MTAB-5061")),                                                           # 51.8%
      S06 = list(
        train = c("E-ANND-1","E-HCAD-10","E-CURD-88","E-CURD-122","E-HCAD-8","E-HCAD-9",
                  "E-MTAB-9543"),
        test  = c("E-HCAD-15","E-ANND-2","E-HCAD-1","E-HCAD-35","E-CURD-119","E-MTAB-10553",
                  "E-HCAD-31","E-GEOD-114530","E-MTAB-8410","E-HCAD-11","E-GEOD-125970",
                  "E-GEOD-84465","E-GEOD-130148","E-ANND-3","E-GEOD-124472","E-MTAB-5061",
                  "E-MTAB-9067","E-GEOD-83139")),                                           # 49.5%
      S07 = list(
        train = c("E-HCAD-15","E-ANND-2","E-HCAD-1","E-HCAD-10","E-CURD-119","E-CURD-88",
                  "E-HCAD-31","E-GEOD-114530","E-MTAB-8410","E-HCAD-8","E-HCAD-11","E-HCAD-9",
                  "E-GEOD-125970","E-GEOD-84465","E-GEOD-130148","E-ANND-3","E-GEOD-124472",
                  "E-MTAB-5061","E-MTAB-9067","E-GEOD-83139"),
        test  = c("E-ANND-1","E-HCAD-35","E-MTAB-10553","E-CURD-122","E-MTAB-9543")),      # 49.9%
      S08 = list(
        train = c("E-ANND-1","E-HCAD-35","E-MTAB-10553","E-MTAB-8410","E-CURD-122",
                  "E-GEOD-125970","E-MTAB-9067"),
        test  = c("E-HCAD-15","E-ANND-2","E-HCAD-1","E-HCAD-10","E-CURD-119","E-CURD-88",
                  "E-HCAD-31","E-GEOD-114530","E-HCAD-8","E-HCAD-11","E-HCAD-9","E-MTAB-9543",
                  "E-GEOD-84465","E-GEOD-130148","E-ANND-3","E-GEOD-124472","E-MTAB-5061",
                  "E-GEOD-83139")),                                                          # 51.4%
      S09 = list(
        train = c("E-ANND-1","E-HCAD-1","E-CURD-119","E-GEOD-83139"),
        test  = c("E-HCAD-15","E-ANND-2","E-HCAD-35","E-HCAD-10","E-CURD-88","E-MTAB-10553",
                  "E-HCAD-31","E-GEOD-114530","E-MTAB-8410","E-CURD-122","E-HCAD-8","E-HCAD-11",
                  "E-HCAD-9","E-MTAB-9543","E-GEOD-125970","E-GEOD-84465","E-GEOD-130148",
                  "E-ANND-3","E-GEOD-124472","E-MTAB-5061","E-MTAB-9067")),                # 49.3%
      S10 = list(
        train = c("E-HCAD-15","E-ANND-2","E-HCAD-35","E-CURD-119","E-CURD-88","E-MTAB-10553",
                  "E-HCAD-31","E-GEOD-114530","E-CURD-122","E-HCAD-8","E-HCAD-11","E-HCAD-9",
                  "E-MTAB-9543","E-GEOD-125970","E-GEOD-84465","E-GEOD-130148","E-ANND-3",
                  "E-GEOD-124472","E-MTAB-5061","E-MTAB-9067"),
        test  = c("E-ANND-1","E-HCAD-1","E-HCAD-10","E-MTAB-8410","E-GEOD-83139"))        # 48.8%
    ),
    "70%" = list(
      S01 = list(
        train = c("E-ANND-1","E-HCAD-15","E-HCAD-10","E-CURD-119","E-CURD-88","E-MTAB-10553",
                  "E-GEOD-114530","E-MTAB-8410","E-GEOD-125970","E-GEOD-84465","E-GEOD-130148",
                  "E-ANND-3","E-GEOD-124472","E-GEOD-83139"),
        test  = c("E-ANND-2","E-HCAD-1","E-HCAD-35","E-HCAD-31","E-CURD-122","E-HCAD-8",
                  "E-HCAD-11","E-HCAD-9","E-MTAB-9543","E-MTAB-5061","E-MTAB-9067")),      # 68.0%
      S02 = list(
        train = c("E-ANND-1","E-ANND-2","E-HCAD-1","E-HCAD-35","E-HCAD-31","E-HCAD-11",
                  "E-HCAD-9","E-MTAB-9543","E-MTAB-5061"),
        test  = c("E-HCAD-15","E-HCAD-10","E-CURD-119","E-CURD-88","E-MTAB-10553",
                  "E-GEOD-114530","E-MTAB-8410","E-CURD-122","E-HCAD-8","E-GEOD-125970",
                  "E-GEOD-84465","E-GEOD-130148","E-ANND-3","E-GEOD-124472","E-MTAB-9067",
                  "E-GEOD-83139")),                                                          # 69.6%
      S03 = list(
        train = c("E-ANND-1","E-HCAD-1","E-HCAD-35","E-HCAD-10","E-CURD-119","E-CURD-88",
                  "E-HCAD-31","E-MTAB-8410","E-CURD-122","E-HCAD-8","E-HCAD-11","E-HCAD-9",
                  "E-MTAB-9543","E-GEOD-125970","E-GEOD-130148","E-ANND-3","E-GEOD-124472",
                  "E-MTAB-5061","E-MTAB-9067"),
        test  = c("E-HCAD-15","E-ANND-2","E-MTAB-10553","E-GEOD-114530","E-GEOD-84465",
                  "E-GEOD-83139")),                                                          # 70.1%
      S04 = list(
        train = c("E-ANND-1","E-HCAD-15","E-ANND-2","E-MTAB-10553","E-CURD-122","E-MTAB-9543",
                  "E-GEOD-84465","E-MTAB-9067"),
        test  = c("E-HCAD-1","E-HCAD-35","E-HCAD-10","E-CURD-119","E-CURD-88","E-HCAD-31",
                  "E-GEOD-114530","E-MTAB-8410","E-HCAD-8","E-HCAD-11","E-HCAD-9",
                  "E-GEOD-125970","E-GEOD-130148","E-ANND-3","E-GEOD-124472","E-MTAB-5061",
                  "E-GEOD-83139")),                                                          # 70.9%
      S05 = list(
        train = c("E-ANND-1","E-HCAD-15","E-HCAD-1","E-HCAD-10","E-MTAB-10553","E-GEOD-114530",
                  "E-HCAD-8","E-HCAD-9","E-GEOD-124472","E-MTAB-5061","E-MTAB-9067",
                  "E-GEOD-83139"),
        test  = c("E-ANND-2","E-HCAD-35","E-CURD-119","E-CURD-88","E-HCAD-31","E-MTAB-8410",
                  "E-CURD-122","E-HCAD-11","E-MTAB-9543","E-GEOD-125970","E-GEOD-84465",
                  "E-GEOD-130148","E-ANND-3")),                                              # 69.7%
      S06 = list(
        train = c("E-ANND-1","E-ANND-2","E-HCAD-35","E-CURD-119","E-CURD-88","E-HCAD-31",
                  "E-MTAB-8410","E-CURD-122","E-HCAD-8","E-HCAD-11","E-HCAD-9","E-GEOD-125970",
                  "E-ANND-3"),
        test  = c("E-HCAD-15","E-HCAD-1","E-HCAD-10","E-MTAB-10553","E-GEOD-114530",
                  "E-MTAB-9543","E-GEOD-84465","E-GEOD-130148","E-GEOD-124472","E-MTAB-5061",
                  "E-MTAB-9067","E-GEOD-83139")),                                            # 71.5%
      S07 = list(
        train = c("E-ANND-1","E-ANND-2","E-HCAD-1","E-HCAD-10","E-CURD-88","E-GEOD-114530",
                  "E-CURD-122","E-HCAD-11","E-MTAB-9543","E-GEOD-125970","E-GEOD-84465",
                  "E-GEOD-130148","E-MTAB-5061","E-MTAB-9067","E-GEOD-83139"),
        test  = c("E-HCAD-15","E-HCAD-35","E-CURD-119","E-MTAB-10553","E-HCAD-31","E-MTAB-8410",
                  "E-HCAD-8","E-HCAD-9","E-ANND-3","E-GEOD-124472")),                      # 70.3%
      S08 = list(
        train = c("E-ANND-1","E-HCAD-15","E-HCAD-35","E-CURD-119","E-MTAB-10553","E-HCAD-31",
                  "E-MTAB-8410","E-CURD-122","E-GEOD-84465","E-GEOD-130148","E-ANND-3",
                  "E-GEOD-124472","E-GEOD-83139"),
        test  = c("E-ANND-2","E-HCAD-1","E-HCAD-10","E-CURD-88","E-GEOD-114530","E-HCAD-8",
                  "E-HCAD-11","E-HCAD-9","E-MTAB-9543","E-GEOD-125970","E-MTAB-5061",
                  "E-MTAB-9067")),                                                            # 70.0%
      S09 = list(
        train = c("E-ANND-1","E-ANND-2","E-HCAD-10","E-CURD-88","E-MTAB-10553","E-HCAD-31",
                  "E-GEOD-114530","E-MTAB-8410","E-HCAD-8","E-HCAD-11","E-MTAB-9543",
                  "E-GEOD-84465","E-GEOD-130148","E-GEOD-124472","E-MTAB-9067"),
        test  = c("E-HCAD-15","E-HCAD-1","E-HCAD-35","E-CURD-119","E-CURD-122","E-HCAD-9",
                  "E-GEOD-125970","E-ANND-3","E-MTAB-5061","E-GEOD-83139")),                # 68.1%
      S10 = list(
        train = c("E-ANND-1","E-HCAD-15","E-HCAD-1","E-HCAD-35","E-HCAD-8","E-HCAD-9",
                  "E-GEOD-125970","E-GEOD-124472","E-MTAB-5061","E-GEOD-83139"),
        test  = c("E-ANND-2","E-HCAD-10","E-CURD-119","E-CURD-88","E-MTAB-10553","E-HCAD-31",
                  "E-GEOD-114530","E-MTAB-8410","E-CURD-122","E-HCAD-11","E-MTAB-9543",
                  "E-GEOD-84465","E-GEOD-130148","E-ANND-3","E-MTAB-9067"))                 # 69.7%
    ),
    "80%" = list(
      S01 = list(
        train = c("E-ANND-1","E-ANND-2","E-HCAD-1","E-HCAD-35","E-HCAD-10","E-CURD-119",
                  "E-CURD-88","E-HCAD-31","E-CURD-122","E-HCAD-8","E-HCAD-9","E-MTAB-9543",
                  "E-GEOD-84465","E-GEOD-130148","E-GEOD-124472","E-MTAB-5061","E-MTAB-9067"),
        test  = c("E-HCAD-15","E-MTAB-10553","E-GEOD-114530","E-MTAB-8410","E-HCAD-11",
                  "E-GEOD-125970","E-ANND-3","E-GEOD-83139")),                              # 79.5%
      S02 = list(
        train = c("E-ANND-1","E-HCAD-15","E-ANND-2","E-HCAD-1","E-MTAB-10553","E-GEOD-114530",
                  "E-MTAB-8410","E-HCAD-11","E-MTAB-9543","E-ANND-3","E-GEOD-83139"),
        test  = c("E-HCAD-35","E-HCAD-10","E-CURD-119","E-CURD-88","E-HCAD-31","E-CURD-122",
                  "E-HCAD-8","E-HCAD-9","E-GEOD-125970","E-GEOD-84465","E-GEOD-130148",
                  "E-GEOD-124472","E-MTAB-5061","E-MTAB-9067")),                            # 79.8%
      S03 = list(
        train = c("E-ANND-1","E-HCAD-15","E-HCAD-35","E-HCAD-10","E-CURD-119","E-CURD-88",
                  "E-HCAD-31","E-GEOD-114530","E-MTAB-8410","E-CURD-122","E-HCAD-8","E-HCAD-11",
                  "E-HCAD-9","E-GEOD-125970","E-GEOD-84465","E-ANND-3","E-GEOD-124472",
                  "E-MTAB-5061"),
        test  = c("E-ANND-2","E-HCAD-1","E-MTAB-10553","E-MTAB-9543","E-GEOD-130148",
                  "E-MTAB-9067","E-GEOD-83139")),                                            # 78.2%
      S04 = list(
        train = c("E-ANND-1","E-HCAD-15","E-ANND-2","E-HCAD-35","E-HCAD-10","E-MTAB-10553",
                  "E-GEOD-114530","E-HCAD-8","E-GEOD-125970","E-GEOD-130148","E-MTAB-9067"),
        test  = c("E-HCAD-1","E-CURD-119","E-CURD-88","E-HCAD-31","E-MTAB-8410","E-CURD-122",
                  "E-HCAD-11","E-HCAD-9","E-MTAB-9543","E-GEOD-84465","E-ANND-3",
                  "E-GEOD-124472","E-MTAB-5061","E-GEOD-83139")),                           # 80.9%
      S05 = list(
        train = c("E-ANND-1","E-HCAD-15","E-ANND-2","E-HCAD-1","E-MTAB-10553","E-MTAB-8410",
                  "E-CURD-122","E-HCAD-9","E-MTAB-9543","E-GEOD-125970","E-GEOD-130148",
                  "E-GEOD-83139"),
        test  = c("E-HCAD-35","E-HCAD-10","E-CURD-119","E-CURD-88","E-HCAD-31","E-GEOD-114530",
                  "E-HCAD-8","E-HCAD-11","E-GEOD-84465","E-ANND-3","E-GEOD-124472",
                  "E-MTAB-5061","E-MTAB-9067")),                                             # 79.9%
      S06 = list(
        train = c("E-ANND-1","E-ANND-2","E-HCAD-1","E-HCAD-35","E-HCAD-10","E-CURD-119",
                  "E-CURD-88","E-HCAD-31","E-CURD-122","E-HCAD-8","E-HCAD-11","E-MTAB-9543",
                  "E-GEOD-84465","E-GEOD-130148","E-ANND-3","E-GEOD-124472","E-MTAB-5061",
                  "E-MTAB-9067","E-GEOD-83139"),
        test  = c("E-HCAD-15","E-MTAB-10553","E-GEOD-114530","E-MTAB-8410","E-HCAD-9",
                  "E-GEOD-125970")),                                                          # 80.1%
      S07 = list(
        train = c("E-ANND-1","E-HCAD-15","E-HCAD-1","E-HCAD-10","E-CURD-119","E-CURD-88",
                  "E-MTAB-10553","E-HCAD-31","E-GEOD-114530","E-MTAB-8410","E-HCAD-11",
                  "E-HCAD-9","E-GEOD-125970","E-GEOD-84465","E-ANND-3","E-GEOD-124472",
                  "E-MTAB-5061","E-MTAB-9067","E-GEOD-83139"),
        test  = c("E-ANND-2","E-HCAD-35","E-CURD-122","E-HCAD-8","E-MTAB-9543",
                  "E-GEOD-130148")),                                                          # 78.1%
      S08 = list(
        train = c("E-ANND-1","E-HCAD-15","E-ANND-2","E-HCAD-1","E-HCAD-35"),
        test  = c("E-HCAD-10","E-CURD-119","E-CURD-88","E-MTAB-10553","E-HCAD-31",
                  "E-GEOD-114530","E-MTAB-8410","E-CURD-122","E-HCAD-8","E-HCAD-11","E-HCAD-9",
                  "E-MTAB-9543","E-GEOD-125970","E-GEOD-84465","E-GEOD-130148","E-ANND-3",
                  "E-GEOD-124472","E-MTAB-5061","E-MTAB-9067","E-GEOD-83139")),             # 78.9%
      S09 = list(
        train = c("E-ANND-1","E-HCAD-15","E-HCAD-35","E-CURD-119","E-CURD-88","E-MTAB-10553",
                  "E-HCAD-31","E-GEOD-114530","E-MTAB-8410","E-CURD-122","E-HCAD-8","E-HCAD-11",
                  "E-HCAD-9","E-MTAB-9543","E-GEOD-125970","E-GEOD-84465","E-GEOD-130148",
                  "E-ANND-3","E-GEOD-124472","E-MTAB-5061","E-GEOD-83139"),
        test  = c("E-ANND-2","E-HCAD-1","E-HCAD-10","E-MTAB-9067")),                      # 78.2%
      S10 = list(
        train = c("E-ANND-1","E-HCAD-15","E-ANND-2","E-HCAD-1","E-CURD-119","E-MTAB-8410",
                  "E-HCAD-8","E-HCAD-11","E-MTAB-9067"),
        test  = c("E-HCAD-35","E-HCAD-10","E-CURD-88","E-MTAB-10553","E-HCAD-31",
                  "E-GEOD-114530","E-CURD-122","E-HCAD-9","E-MTAB-9543","E-GEOD-125970",
                  "E-GEOD-84465","E-GEOD-130148","E-ANND-3","E-GEOD-124472","E-MTAB-5061",
                  "E-GEOD-83139"))                                                            # 79.3%
    )
  )
)
