"""Train and audit a shared spectral conditional generator on synthetic fields.

Training uses noisy field coefficients, not their conditional means. The
independent audit can evaluate Gaussian conditional W2 because this is a known
synthetic data-generating process. Statistical certificates are explicitly
separate from empirical coupling costs. CPU float64, fixed seeds, no downloads.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import time
import numpy as np
import scipy
from scipy.special import softmax, zeta
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SEEDS = [20260912, 20260913, 20260914]
TRAIN_N, TRAIN_M, AUDIT_N, AUDIT_M = 8192, 32, 65536, 256
STEPS, BATCH, WIDTH = 1800, 2048, 12
DELTA = .01


def fields(rng, n, m):
    noise = rng.standard_normal((n, m))
    u = np.empty_like(noise)
    u[:, 0] = noise[:, 0]
    for k in range(1, m):
        u[:, k] = (.6*np.tanh(u[:, 0])+.4*np.tanh(u[:, k-1])+noise[:, k])/(k+1)
    return u, noise


class MeanNet(torch.nn.Module):
    def __init__(self, seed, ablate=False):
        super().__init__()
        torch.manual_seed(seed)
        self.rows = torch.nn.Parameter(torch.randn(WIDTH, 2)*.5)
        self.logits = torch.nn.Parameter(torch.zeros(WIDTH))
        self.log_sigma = torch.nn.Parameter(torch.tensor(0.))
        self.ablate = ablate

    def forward(self, x):
        if self.ablate:
            x = x*torch.tensor([1., 0.])
        return torch.tanh(x @ self.rows.T) @ self.logits.softmax(0)

    def constrain(self):
        with torch.no_grad():
            norms = self.rows.norm(dim=1, keepdim=True)
            self.rows.mul_(torch.clamp(1.3/norms, max=1.))
            self.log_sigma.clamp_(np.log(.7), np.log(1.3))

    def record(self):
        a = self.rows.detach().numpy().copy()
        if self.ablate:
            a[:, 1] = 0
        w = self.logits.detach().softmax(0).numpy()
        return {"rows": a.tolist(), "weights": w.tolist(),
                "sigma": float(self.log_sigma.detach().exp()),
                "global_mean_bound": 1.,
                "global_context_lipschitz": float(w @ np.linalg.norm(a, axis=1)),
                "parameters": sum(p.numel() for p in self.parameters()),
                "ablation": self.ablate}


def mean_numpy(x, record):
    a, w = np.asarray(record["rows"]), np.asarray(record["weights"])
    return (np.tanh(x[:, 0, None]*a[None, :, 0]+x[:, 1, None]*a[None, :, 1])*w).sum(axis=1)


def train(seed, ablate=False):
    rng = np.random.default_rng(seed)
    u, _ = fields(rng, TRAIN_N, TRAIN_M)
    model = MeanNet(seed, ablate)
    opt = torch.optim.Adam(model.parameters(), lr=.012)
    loss_log = []
    for step in range(STEPS):
        i = rng.integers(TRAIN_N, size=BATCH)
        # Half the pairs use the first refinement, the rest are uniform later bands.
        j = rng.integers(2, TRAIN_M, size=BATCH)
        j[rng.random(BATCH)<.5] = 1
        x = torch.from_numpy(np.stack([u[i, 0], u[i, j-1]], axis=1))
        y = torch.from_numpy((j+1)*u[i, j])
        resid = y-model(x)
        loss = .5*(resid.square()*torch.exp(-2*model.log_sigma)).mean()+model.log_sigma
        opt.zero_grad()
        loss.backward()
        opt.step()
        model.constrain()
        if step % 100 == 0 or step == STEPS-1:
            loss_log.append([step, float(loss.detach())])
    return model.record(), loss_log


def generate(noise, record):
    out = np.empty_like(noise)
    out[:, 0] = noise[:, 0]
    for k in range(1, noise.shape[1]):
        x = np.stack([out[:, 0], out[:, k-1]], axis=1)
        out[:, k] = (mean_numpy(x, record)+record["sigma"]*noise[:, k])/(k+1)
    return out


def audit(record, target, noise, multiplicity):
    n, m = target.shape
    conditional = np.zeros(m)
    upper = np.zeros(m)
    log_term = np.log(2*multiplicity*(m-1)/DELTA)
    for k in range(1, m):
        x = np.stack([target[:, 0], target[:, k-1]], axis=1)
        truth = .6*np.tanh(x[:, 0])+.4*np.tanh(x[:, 1])
        error = (truth-mean_numpy(x, record))**2
        assert np.min(error) >= 0 and np.max(error) <= 4+1e-12
        conditional[k] = (error.mean()+(1-record["sigma"])**2)/(k+1)**2
        # Maurer--Pontil (2009), Theorem 4; sample variance uses ddof=1.
        bonus = np.sqrt(2*error.var(ddof=1)*log_term/n)+7*4*log_term/(3*(n-1))
        upper[k] = (min(4., error.mean()+bonus)+(1-record["sigma"])**2)/(k+1)**2
    L = record["global_context_lipschitz"]/np.arange(1, m+1)
    L[0] = 0
    L[1] *= np.sqrt(2)  # Context extraction duplicates x_1 at j=2.
    d_emp, d_high = np.zeros(m), np.zeros(m)
    for k in range(1, m):
        d_emp[k] = np.hypot(d_emp[k-1], np.sqrt(conditional[k])+L[k]*d_emp[k-1])
        d_high[k] = np.hypot(d_high[k-1], np.sqrt(upper[k])+L[k]*d_high[k-1])
    generated = generate(noise, record)
    rows = []
    for retained in [8, 16, 32, 64, 128, 256]:
        projected_squared = ((target[:, :retained]-generated[:, :retained])**2).sum(axis=1)
        full_squared = projected_squared+(target[:, retained:]**2).sum(axis=1)
        tail = 2*zeta(2, retained+1)
        rows.append({"m": retained,
            "projected_coupling_rms": float(np.sqrt(projected_squared.mean())),
            "finite_reference_coupling_rms": float(np.sqrt(full_squared.mean())),
            "squared_coupling_standard_error": float(full_squared.std(ddof=1)/np.sqrt(n)),
            "estimated_recurrence": float(d_emp[retained-1]),
            "confidence_recurrence": float(d_high[retained-1]),
            "confidence_full_law_bound": float(np.sqrt(tail+d_high[retained-1]**2)),
            "target_tail_bound": float(np.sqrt(tail))})
    prefix = generate(noise[:64, :32], record)
    assert np.array_equal(prefix, generated[:64, :32])
    assert np.all(np.sqrt(conditional) <= np.sqrt(upper)+1e-14)
    # These empirical inequalities are diagnostics, not the confidence theorem.
    return {"rows": rows, "conditional_w2_squared_estimates": conditional.tolist(),
            "conditional_w2_squared_upper": upper.tolist(), "L": L.tolist(),
            "prefix_max_difference": float(np.max(abs(prefix-generated[:64, :32]))),
            "aggregate_conditional_error": float(np.sqrt(conditional.sum()))}


def main():
    torch.set_default_dtype(torch.float64)
    torch.set_num_threads(1)
    torch.use_deterministic_algorithms(True)
    for d in ["results", "figures", "tables", "models"]:
        (ROOT/d).mkdir(exist_ok=True)
    started = time.monotonic()
    fitted = []
    for seed in SEEDS:
        record, log = train(seed)
        fitted.append({"seed": seed, "model": record, "training_loss": log})
        print(f"trained seed {seed}: L={record['global_context_lipschitz']:.6f}, sigma={record['sigma']:.6f}", flush=True)
    ablation, ablation_log = train(SEEDS[0], ablate=True)
    target, noise = fields(np.random.default_rng(937105), AUDIT_N, AUDIT_M)
    for fit in fitted:
        fit["audit"] = audit(fit["model"], target, noise, len(SEEDS))
    ablation_audit = audit(ablation, target, noise, len(SEEDS)+1)
    # Mean-zero independent Gaussian baseline, with the correct innovation scale.
    independent = {"rows": [[0., 0.]], "weights": [1.], "sigma": 1.,
                   "global_context_lipschitz": 0., "global_mean_bound": 0.}
    independent_out = generate(noise, independent)
    baseline_rms = np.sqrt(((target-independent_out)**2).sum(axis=1).mean())
    report = {"kind": "Trained synthetic spectral conditional neural generator",
        "training": {"independent_fields": TRAIN_N, "retained_modes": TRAIN_M,
            "steps": STEPS, "batch_pairs": BATCH, "width": WIDTH,
            "objective": "Gaussian conditional negative log-likelihood, constant removed",
            "context_sampling": "half mode 2; half uniform modes 3 through 32",
            "optimizer": "Adam, learning rate 0.012, row norms projected to at most 1.3"},
        "validation": {"independent_fields": AUDIT_N, "retained_modes": AUDIT_M,
            "seed": 937105, "confidence": 1-DELTA,
            "scope": "simultaneously all 3 full-context fitted networks and modes 2..256",
            "method": "bounded empirical Bernstein, union bound; synthetic conditional means are known",
            "held_out_data_used_for_selection": False},
        "fits": fitted,
        "ablation": {"description": "remove previous-band input; retrain on the same noisy field pairs",
            "model": ablation, "training_loss": ablation_log, "audit": ablation_audit},
        "independent_gaussian_baseline_projected_rms": float(baseline_rms),
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
            "torch": torch.__version__, "scipy": scipy.__version__, "dtype": "float64", "torch_threads": 1}}
    (ROOT/"results/learned_operator.json").write_text(json.dumps(report, indent=2)+"\n")
    for fit in fitted:
        (ROOT/f"models/spectral_generator_{fit['seed']}.json").write_text(json.dumps(fit["model"],indent=2)+"\n")
    (ROOT/"models/spectral_generator_ablation.json").write_text(json.dumps(ablation,indent=2)+"\n")
    lines = [r"\begin{tabular}{rrrrr}", r"\toprule",
        r"Modes & Coupling RMS & Estimated $d_m$ & 99\% $d_m$ & 99\% full-law bound \\", r"\midrule"]
    for row in fitted[0]["audit"]["rows"]:
        lines.append(f"{row['m']} & {row['finite_reference_coupling_rms']:.5f} & {row['estimated_recurrence']:.5f} & {row['confidence_recurrence']:.5f} & {row['confidence_full_law_bound']:.5f} "+r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (ROOT/"tables/learned_operator.tex").write_text("\n".join(lines)+"\n")
    plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False})
    fig, axes = plt.subplots(1,2,figsize=(10,3.7),layout="constrained")
    for k,fit in enumerate(fitted):
        rows = fit["audit"]["rows"]
        axes[0].loglog([r["m"] for r in rows],[r["finite_reference_coupling_rms"] for r in rows],"o-",alpha=.7,label=f"Learned, seed {k+1}")
    rows=fitted[0]["audit"]["rows"]
    axes[0].loglog([r["m"] for r in rows],[r["confidence_full_law_bound"] for r in rows],"k--",label="99% bound, seed 1")
    axes[0].axvline(TRAIN_M,color="grey",ls=":")
    axes[0].set(xlabel="Retained modes (trained through 32)",ylabel="Coupling RMS / law bound",title="Shared network across resolutions")
    axes[0].legend(fontsize=8)
    vals=[fit["audit"]["rows"][-1]["projected_coupling_rms"] for fit in fitted]
    vals += [ablation_audit["rows"][-1]["projected_coupling_rms"],float(baseline_rms)]
    axes[1].bar(range(5),vals,color=["#23658a"]*3+["#bd583b","#596574"])
    axes[1].set_xticks(range(5),["Seed 1","Seed 2","Seed 3","No last\nband","No\ncontext"])
    axes[1].set(ylabel="256-mode coupling RMS",title="Context ablations")
    fig.savefig(ROOT/"figures/learned_operator.png",dpi=220,bbox_inches="tight",facecolor="white")
    plt.close(fig)
    print(json.dumps({"seconds":time.monotonic()-started,
        "learned_rms": vals[:3], "ablated_rms": vals[3:],
        "first_seed_bound": rows[-1]["confidence_full_law_bound"]}),flush=True)


if __name__ == "__main__":
    main()
