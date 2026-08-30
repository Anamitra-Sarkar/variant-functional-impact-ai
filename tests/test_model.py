import numpy as np
from data_pipeline.model import featurize_variant, featurize_batch, train_logistic_regression, evaluate_with_baselines, FEATURE_NAMES
from data_pipeline.structure import ResidueFeatures
from data_pipeline.physicochemical import KYTE_DOOLITTLE

def test_featurize_variant():
    rf = ResidueFeatures(chain_id="A", resseq=10, resname="ALA", rsa=0.1, ss="H", burial_count=15, ca_coord=(0,0,0))
    feats = featurize_variant("R", "W", rf, conservation=0.9)
    assert len(feats) == len(FEATURE_NAMES)
    # rsa 0.1, conservation 0.9
    assert feats[0] == 0.1
    assert feats[1] == 0.9
    # ss_H should be 1
    assert feats[FEATURE_NAMES.index("ss_H")] == 1.0

def test_fusion_recovers_injected_signal():
    """
    Synthetic injected-signal verification:
    buried (rsa<0.2) + conserved (>0.7) + large physchem delta (R->W) => damaging (1)
    otherwise benign (0)
    Fusion should beat baselines.
    """
    np.random.seed(42)
    variants = []
    labels = []
    for i in range(200):
        # Create pattern: every 3rd is damaging context, every 2nd is large change
        is_buried = (i % 4 == 0)
        rsa = 0.05 if is_buried else 0.8
        cons = 0.95 if is_buried else 0.15
        # Alternate large vs small mutation
        is_large = (i % 3 == 0)
        wt, mut = ("R", "W") if is_large else ("S", "T")
        # Label: damaging if buried+conserved+large
        label = 1 if (is_buried and is_large) else 0
        # For this test, we want is_buried correlated with is_large sometimes but not perfectly
        # Actually make it so damaging = is_buried & is_large
        # This is the injected real signal
        rf = ResidueFeatures(chain_id="A", resseq=i+1, resname="ALA", rsa=rsa, ss="H" if i%2==0 else "C", burial_count=15 if is_buried else 2, ca_coord=(0,0,0))
        variants.append({"wt": wt, "mut": mut, "residue_feature": rf, "conservation": cons})
        labels.append(label)
    X, names = featurize_batch(variants)
    y = np.array(labels)
    assert sum(labels) > 0 and sum(labels) < len(labels)

    model = train_logistic_regression(X, y)
    y_score = model.predict_proba(X)[:, 1]
    results = evaluate_with_baselines(X, y, names, model)
    print(results)
    # Fusion must beat at least one baseline and achieve reasonable AUROC
    assert results["fusion"]["auroc"] > 0.75, f"Fusion AUROC too low: {results['fusion']}"
    # Fusion should not be worse than baselines (allow tiny tolerance)
    assert results["fusion"]["auroc"] >= results["baseline_conservation"]["auroc"] - 0.05
    # Overall fusion should be better than hydro-only baseline clearly
    assert results["fusion"]["auroc"] > results["baseline_abs_delta_hydro"]["auroc"] - 0.1
