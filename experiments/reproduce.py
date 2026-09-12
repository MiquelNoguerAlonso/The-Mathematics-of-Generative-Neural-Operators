#!/usr/bin/env python3
"""Reproduce the paper's analytic and synthetic numerical illustrations.

These are oracle-velocity / explicit-kernel experiments, not trained-network
benchmarks. Every reported Wasserstein value is either exact for Gaussians or
explicitly labelled as a bound or a Monte Carlo coupling estimate.
"""
from __future__ import annotations

import json
from pathlib import Path
import platform
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.special import polygamma
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
TAB = ROOT / "tables"
RES = ROOT / "results"
for directory in (FIG, TAB, RES):
    directory.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 12,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelsize": 11.5, "axes.titlesize": 12,
    "legend.fontsize": 10, "figure.dpi": 130,
    "savefig.dpi": 230,
})
BLUE, RED, GREEN, GOLD = "#235789", "#b64545", "#24836b", "#c28b28"


def psd_sqrt(a):
    a = (a + a.T) / 2
    vals, vecs = np.linalg.eigh(a)
    if vals.min() < -1e-8:
        raise ValueError("Matrix is not positive semidefinite")
    return (vecs * np.sqrt(np.maximum(vals, 0))) @ vecs.T


def gaussian_w2(a, b):
    rb = psd_sqrt(b)
    squared = np.trace(a + b - 2 * psd_sqrt(rb @ a @ rb))
    if abs(squared) < 2e-13:
        squared = 0.0
    return float(np.sqrt(max(squared, 0.0)))


def full_velocity(s, rho):
    sigma = np.array([[1.0, rho], [rho, 1.0]])
    cov = (1-s)**2 * np.eye(2) + s*s * sigma
    return (s * sigma - (1-s)*np.eye(2)) @ np.linalg.inv(cov)


def triangular_velocity(s, rho):
    a = full_velocity(s, rho)
    d = (1-s)**2 + s*s
    a[0, :] = [(2*s-1)/d, 0]
    return a


def defect(s, rho):
    d = (1-s)**2 + s*s
    return rho*rho*s*s*(1-s)**2 / (d*(d*d-s**4*rho*rho))


def terminal_map(rho, triangular=False):
    field = triangular_velocity if triangular else full_velocity
    def rhs(s, state):
        return (field(s, rho) @ state.reshape(2, 2)).ravel()
    sol = solve_ivp(rhs, (0, 1), np.eye(2).ravel(), method="DOP853",
                    rtol=1e-11, atol=1e-13)
    assert sol.success
    return sol.y[:, -1].reshape(2, 2)


def gaussian_experiments():
    rows = []
    s_grid = np.linspace(0, 1, 401)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), constrained_layout=True)
    for rho, color in zip((.25, .5, .8, .95), (BLUE, GREEN, GOLD, RED)):
        sigma = np.array([[1., rho], [rho, 1.]])
        f_full = terminal_map(rho)
        f_tri = terminal_map(rho, triangular=True)
        f_ref = np.array([[1., 0.], [rho, np.sqrt(1-rho*rho)]])
        cov_tri = f_tri @ f_tri.T
        integral, numerical_err = quad(lambda s: defect(s, rho), 0, 1,
                                       epsabs=2e-12, epsrel=2e-12)
        assert np.max(np.abs(f_full@f_full.T-sigma)) < 5e-10
        assert np.max(np.abs(f_ref@f_ref.T-sigma)) < 2e-15
        assert abs(defect(.5, rho)-2*rho*rho/(4-rho*rho)) < 2e-15
        for s in (.1, .25, .5, .8, .95):
            cov = (1-s)**2*np.eye(2)+s*s*sigma
            condvar = cov[1,1]-cov[1,0]**2/cov[0,0]
            assert abs(full_velocity(s, rho)[0,1]**2*condvar-defect(s, rho)) < 2e-12
        row = dict(rho=rho, defect_midpoint=defect(.5,rho),
                   integrated_defect=integral, quadrature_error=numerical_err,
                   triangular_w2=gaussian_w2(sigma,cov_tri),
                   triangular_covariance=cov_tri.tolist(),
                   triangular_correlation=float(cov_tri[0,1]/np.sqrt(cov_tri[0,0]*cov_tri[1,1])),
                   full_covariance_residual=float(np.max(np.abs(f_full@f_full.T-sigma))),
                   triangular_transport_cost=float(2-2*np.sqrt(1-rho*rho)),
                   unrestricted_transport_cost=float(4-2*(np.sqrt(1+rho)+np.sqrt(1-rho))))
        row['transport_premium']=row['triangular_transport_cost']-row['unrestricted_transport_cost']
        rows.append(row)
        axes[0].plot(s_grid, defect(s_grid,rho), color=color, label=rf"$\rho={rho}$")
    axes[0].set(xlabel="Transport time $s$", ylabel="Irreducible squared velocity error",
                title="(a) Irreducible velocity defect")
    axes[0].legend(frameon=False)
    axes[1].plot([r['rho'] for r in rows], [r['triangular_w2'] for r in rows],
                 'o-', color=RED, label="Projected global velocity")
    axes[1].plot([r['rho'] for r in rows], np.zeros(len(rows)), 's--',
                 color=GREEN, label="Conditional refinement (exact)")
    axes[1].set(xlabel=r"Target correlation $\rho$", ylabel="Endpoint $W_2$ error",
                title="(b) Endpoint distribution error")
    axes[1].legend(frameon=False)
    fig.savefig(FIG / "bridge_obstruction.png")
    plt.close(fig)
    lines = [r"\begin{tabular}{rrrrr}", r"\toprule",
             r"$\rho$ & $D_P(1/2)$ & $\int_0^1D_P(s)\,ds$ & Endpoint $W_2$ & Correlation \\",
             r"\midrule"]
    for r in rows:
        lines.append(f"{r['rho']:.2f} & {r['defect_midpoint']:.6f} & {r['integrated_defect']:.6f} & {r['triangular_w2']:.6f} & {r['triangular_correlation']:.6f} " + r"\\")
    lines.extend([r"\bottomrule",r"\end{tabular}"])
    (TAB / "gaussian.tex").write_text("\n".join(lines)+"\n")

    rho = .8
    sigma = np.array([[1.,rho],[rho,1.]])
    tri = terminal_map(rho, True)
    angle = np.linspace(0,2*np.pi,401)
    circle = np.array([np.cos(angle),np.sin(angle)])
    fig, axes = plt.subplots(1,2,figsize=(9,3.6),constrained_layout=True)
    curves = [(psd_sqrt(sigma),BLUE,"Target law"),(tri,RED,"Projected global flow")]
    for f,color,label in curves:
        ellipse = f@circle
        axes[0].plot(*ellipse,color=color,label=label,lw=2)
    axes[0].set(xlabel="Coarse coordinate",ylabel="Detail coordinate",
                title="(a) Covariance contours",aspect="equal")
    axes[0].legend(frameon=False)
    sig = np.sqrt(1-rho*rho)
    for seed, color in zip((-1.,0.,1.),(BLUE,GOLD,RED)):
        path = s_grid*rho + np.sqrt((1-s_grid)**2+s_grid*s_grid*sig*sig)*seed
        axes[1].plot(s_grid,path,color=color,label=f"Initial detail {seed:g}")
    axes[1].set(xlabel="Refinement time $r$",ylabel="Detail coordinate",
                title="(b) Refinement at fixed coarse value")
    axes[1].legend(frameon=False)
    fig.savefig(FIG / "conditional_refinement.png")
    plt.close(fig)
    return rows


def refinement_experiments():
    seed, n, maximum = 20260911, 60000, 256
    rng = np.random.default_rng(seed)
    modes = [2,4,8,16,32,64,128,256]
    b, delta = .08, .12
    first = rng.normal(size=n)
    prev_true, prev_hat = first.copy(), first.copy()
    energy = first*first
    sqerr = np.zeros(n)
    recorded = {}
    certificate = 0.
    samples = [first[:4].copy()]
    context_sum = 0.
    epsilon_sum = 0.
    d_sequence = [0.]
    for j in range(2,maximum+1):
        noise = rng.normal(size=n)
        true = (.6*np.tanh(first)+.4*np.tanh(prev_true)+noise)/j
        generated = (.6*np.tanh(first)+.4*np.tanh(prev_hat)+b+(1+delta)*noise)/j
        energy += true*true
        sqerr += (true-generated)**2
        eps = np.sqrt(b*b+delta*delta)/j
        lip = (1. if j==2 else np.sqrt(.52))/j
        certificate = np.sqrt(certificate**2 + (eps+lip*certificate)**2)
        context_sum += lip*lip
        epsilon_sum += eps*eps
        d_sequence.append(certificate)
        if j in modes:
            recorded[j] = (sqerr.copy(),energy.copy(),certificate)
        if j <=128:
            samples.append(true[:4].copy())
        prev_true, prev_hat = true,generated
    rows = []
    for m in modes:
        partialerr,partialenergy,d = recorded[m]
        cost = partialerr+energy-partialenergy
        mean = float(cost.mean())
        se = float(cost.std(ddof=1)/np.sqrt(n))
        rms = np.sqrt(mean)
        tau2_bound = float(2*polygamma(1,m+1))
        bound = np.sqrt(d*d+tau2_bound)
        assert mean < bound*bound + 5*se
        rows.append(dict(m=m, coupling_rms_to_256_mode_reference=float(rms),
                         coupling_squared_se=se, refinement_certificate=float(d),
                         tail_energy_upper_bound=tau2_bound,
                         full_law_w2_upper_bound=float(bound)))

    fig,axes=plt.subplots(1,2,figsize=(9,3.8),constrained_layout=True)
    axes[0].loglog(modes,[r['full_law_w2_upper_bound'] for r in rows],
                   'o-',color=BLUE,label="Infinite-law certificate")
    axes[0].loglog(modes,[r['coupling_rms_to_256_mode_reference'] for r in rows],
                   's--',color=GREEN,label="256-mode coupling RMS estimate")
    axes[0].loglog(modes,[r['refinement_certificate'] for r in rows],
                   ':',color=RED,label="Refinement contribution")
    axes[0].set(xlabel="Retained modes $m$",ylabel="Error / upper bound",
                title="(a) Nonlinear field certificate")
    axes[0].legend(frameon=False,fontsize=9)
    ms=np.arange(1,257)
    for power,color in zip((1.,1.5,2.),(BLUE,GREEN,GOLD)):
        if power==1.:
            tail=polygamma(1,ms+1)
        else:
            from scipy.special import zeta
            tail=zeta(2*power,ms+1)
        axes[1].loglog(ms,np.sqrt(tail),color=color,label=rf"$\lambda_j=j^{{-{2*power:g}}}$")
    axes[1].set(xlabel="Retained modes $m$",ylabel="Exact truncation $W_2$",
                title="(b) Exact Gaussian truncation error")
    axes[1].legend(frameon=False)
    fig.savefig(FIG / "refinement_certificate.png")
    plt.close(fig)

    coefficients=np.stack(samples,axis=1)
    x=np.linspace(0,1,601)
    basis=np.sqrt(2)*np.sin(np.pi*np.arange(1,129)[:,None]*x[None,:])
    fig,axes=plt.subplots(2,2,figsize=(9,5.6),constrained_layout=True)
    for idx,ax in enumerate(axes.flat):
        for m,color,lw in [(8,GOLD,1.1),(32,GREEN,1.1),(128,BLUE,1.3)]:
            ax.plot(x,coefficients[idx,:m]@basis[:m],color=color,lw=lw,
                    label=f"{m} modes",alpha=.85)
        ax.set(xlabel="Spatial coordinate",ylabel="Field value",title=f"Realization {idx+1}: shared coefficients")
        if idx==0:
            ax.legend(frameon=False,ncol=3)
    fig.savefig(FIG / "coherent_fields.png")
    plt.close(fig)
    lines=[r"\begin{tabular}{rrrr}",r"\toprule",
           r"Modes & Coupling RMS$^{\dagger}$ & $d_m$ & Full-law bound \\",r"\midrule"]
    for r in rows:
        lines.append(f"{r['m']} & {r['coupling_rms_to_256_mode_reference']:.6f} & {r['refinement_certificate']:.6f} & {r['full_law_w2_upper_bound']:.6f} " + r"\\")
    lines.extend([r"\bottomrule",r"\end{tabular}"])
    (TAB / "refinement.tex").write_text("\n".join(lines)+"\n")
    return dict(seed=seed,samples=n,reference_modes=maximum,bias=b,scale_error=delta,
                sum_context_lipschitz_squared=context_sum,
                sum_conditional_errors_squared=epsilon_sum,rows=rows)


def verification():
    checks=[]
    rng=np.random.default_rng(1983)
    # Exact projection Pythagoras with correlated Gaussian coordinates.
    a=rng.normal(size=(6,6)); cov=a@a.T+.2*np.eye(6)
    b=rng.normal(size=(3,3)); small=b@b.T+.3*np.eye(3)
    embed=np.zeros((6,6)); embed[:3,:3]=small
    lhs=gaussian_w2(cov,embed)**2
    rhs=np.trace(cov[3:,3:])+gaussian_w2(cov[:3,:3],small)**2
    assert abs(lhs-rhs)<2e-8
    checks.append({"name":"orthogonal_projection_transport_identity","absolute_residual":abs(lhs-rhs)})
    # Sharp recurrence: deterministic target at zero, radial learned kernels.
    state=[]; d=0.
    for j in range(40):
        e=.03/(j+1); lip=.2/(j+1)**.75
        detail=e+lip*np.linalg.norm(state)
        state.append(detail)
        d=np.sqrt(d*d+(e+lip*d)**2)
    residual=abs(np.linalg.norm(state)-d)
    assert residual<2e-15
    checks.append({"name":"sharp_multilevel_recurrence","absolute_residual":residual})
    # The exact Gaussian derivative and cross-scale commutator criterion.
    for rho in (0.,.25,.8,.95):
        for s in (.1,.37,.5,.91):
            d=(1-s)**2+s*s
            cross=rho*s*(1-s)/(d*d-s**4*rho*rho)
            assert abs(full_velocity(s,rho)[0,1]-cross)<2e-13
    checks.append({"name":"gaussian_obstruction_algebra","status":"passed"})
    # Conditional Gaussian refinement solves its own population ODE.
    rho=.8; sig2=1-rho*rho; x=1.3; z=-.7
    def rhs(r,y):
        a=(1-r)**2+r*r*sig2
        return [rho*x+(r*sig2-(1-r))/a*(y[0]-r*rho*x)]
    sol=solve_ivp(rhs,(0,1),[z],method="DOP853",rtol=1e-11,atol=1e-13)
    residual=abs(sol.y[0,-1]-(rho*x+np.sqrt(sig2)*z))
    assert residual<5e-10
    checks.append({"name":"conditional_flow_endpoint","absolute_residual":float(residual)})
    # Optimal triangular transport and its kinetic-action premium, in six dimensions.
    a=rng.normal(size=(6,6)); sigma=a@a.T+.5*np.eye(6)
    lower=np.linalg.cholesky(sigma)
    j_tri=float(np.linalg.norm(lower-np.eye(6),'fro')**2)
    j_formula=float(np.trace(sigma)+6-2*np.trace(lower))
    j_opt=gaussian_w2(np.eye(6),sigma)**2
    positive_gap=float(np.trace((lower-psd_sqrt(sigma)).T @
                               np.linalg.inv(psd_sqrt(sigma)) @ (lower-psd_sqrt(sigma))))
    assert abs(j_tri-j_formula)<5e-13
    assert abs((j_tri-j_opt)-positive_gap)<5e-12
    for s in (.1,.5,.9):
        velocity=(lower-np.eye(6))@np.linalg.inv((1-s)*np.eye(6)+s*lower)
        assert np.max(np.abs(np.triu(velocity,1)))<1e-13
    checks.append({"name":"triangular_transport_premium_and_autonomous_bridge",
                   "absolute_residual":abs((j_tri-j_opt)-positive_gap)})
    return checks


def main():
    checks=verification()
    gaussian=gaussian_experiments()
    refinement=refinement_experiments()
    result=dict(environment={"python":platform.python_version(),"numpy":np.__version__,
                             "scipy":scipy.__version__,"matplotlib":matplotlib.__version__},
                verification=checks,gaussian=gaussian,refinement=refinement)
    (RES/"metrics.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"verification":checks,"gaussian":gaussian,
                      "refinement":refinement},indent=2))


if __name__=="__main__":
    main()
