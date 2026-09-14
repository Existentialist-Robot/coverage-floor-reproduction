"""
Substrate + instrument for the COMPLEX NETWORKS 2026 simulation.

Implements the frozen spec (sim-spec-FROZEN-2026-08-04.md) sections 2-4:

  SUBSTRATE (deliberately standard):
    N participant agents; programs arrive in cycles; each program executes the
    six-friction-point lifecycle of main-draft-v1.md section 4.2 --
      1 challenge intake        -> node arrival + attribute assignment
      2 team formation          -> edge / motif formation (Guimera newcomer-vs-repeat-tie)
      3 sprint progression      -> intra-team edge dynamics (stalls, iterations)
      4 evaluation              -> edge weighting (evaluator concordance)
      5 escalation              -> inter-layer edge creation (time-to-decision)
      6 partner re-engagement   -> repeat edge / temporal persistence
    Output: temporal multilayer edge log (team / escalation / partner layers).
    Latent adaptive capacity A(t) drives repeat-tie formation propensity; the
    windowed cumulative repeat-tie graph therefore percolates, and the giant
    component collapses as A(t) declines (degradation scenarios).  Control
    scenarios hold A constant.

  INSTRUMENT (the novel layer) -- four dials:
    1 consent coverage c, arms: uniform (MCAR, reference) and
      centrality-correlated (MNAR, primary); revocation makes it time-varying MNAR.
    2 proxy error eps: white / red (AR(1)) / systematic, three severities each.
    3 revocation r: retroactive removal of already-recorded observations,
      propagating through published aggregates (the observer re-derives its whole
      history from the currently-consented set).
    4 reward coupling g (coverage-and-coupling campaign only): agents respond to the RECORDED proxy under
      partial coverage, never to the true state; includes a Bol-style
      withdrawal/attrition channel.

Pure Python stdlib + numpy.  No LLM calls, no network access.  Every run is
keyed by an explicit integer seed and is bit-for-bit reproducible.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, asdict
from typing import Any

import numpy as np

# --------------------------------------------------------------------------
# RNG stream separation.
#
# Each subsystem draws from its own independent stream spawned from one
# SeedSequence.  This is what makes "matched seeds" meaningful across arms:
# switching an instrument dial off silences that dial's stream without
# re-phasing the substrate streams.
# --------------------------------------------------------------------------
STREAMS = (
    "intake",      # node arrival, challenge attributes
    "formation",   # Guimera team assembly
    "progress",    # stalls, iterations, evaluation concordance
    "escalate",    # escalation decisions
    "reengage",    # repeat-tie formation / partner return
    "consent",     # consent assignment + revocation draws
    "proxy",       # measurement error
    "coupling",    # reward-coupling response + attrition (coverage-and-coupling campaign)
)


def make_streams(seed: int) -> dict[str, np.random.Generator]:
    ss = np.random.SeedSequence(int(seed))
    children = ss.spawn(len(STREAMS))
    return {name: np.random.Generator(np.random.PCG64(child))
            for name, child in zip(STREAMS, children)}


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
@dataclass
class Config:
    # ---- substrate ----
    n_agents: int = 400
    n_cycles: int = 320                # RECORDED cycles
    burn_in: int = 170                 # unrecorded cycles, so the recorded series
                                       # starts on the equilibrated branch
    warm_start: int = 40               # first cycles of burn-in bypass the structural
                                       # feedback term (p_re = A) so the graph does not
                                       # have to escape the empty state by itself;
                                       # a warm start, not part of the studied dynamics
    programs_per_cycle: int = 12
    team_size_choices: tuple[int, ...] = (3, 4, 5)
    program_duration: int = 6          # cycles from open to close (lagging metric lag)
    tie_window: int = 20               # cycles a repeat tie stays in the cumulative graph
                                       # (also the structural correlation time; kept
                                       # shorter than the EWS detrending window so
                                       # in-window detrending does not eat the signal)
    p_newcomer: float = 0.28           # Guimera p: newcomer fraction of team slots
    p_repeat_tie: float = 0.55         # Guimera q: incumbent slot filled from a prior collaborator
    turnover: bool = True              # recycle the longest-inactive agent into a fresh
                                       # newcomer when the newcomer reservoir is empty,
                                       # so the Guimera p channel stays live all run

    # ---- latent adaptive capacity A(t) ----
    # Re-engagement propensity is A(t) * g(z), with z the current mean degree of
    # the windowed cumulative repeat-tie graph and g a cooperative (Hill) term:
    # coordination begets coordination.  The resulting mean-field map z = K(A) g(z)
    # has a FOLD (saddle-node) bifurcation, so the collapse of the giant component
    # is bifurcation-like by construction and exhibits critical slowing down.
    # Bifurcation validation measures the fold location, hysteresis and finite-size behaviour.
    reengage_theta: float = 0.05       # baseline re-engagement with no structure
    reengage_hill: float = 1.2         # h: half-saturation mean degree
    a_high: float = 0.48
    a_low: float = 0.31
    degrade_start: int = 50            # A(t) constant before this recorded cycle
    scenario: str = "degrade"          # "degrade" | "control"
    a_control: float = 0.48            # A held here in control scenario
    burn_in_a: float | None = None     # A during burn-in; None => a_high (degrade)
                                       # / a_control (control).  Setting it low vs
                                       # high is how Bifurcation validation probes hysteresis.

    # ---- observables ----
    formation_base: float = 4.0        # L0
    formation_kappa: float = 9.0       # kappa (search cost, divided by 1 + s_norm)
    formation_sigma: float = 0.35      # idiosyncratic per-program noise (true process)
    escalation_base: float = 3.0       # E0
    escalation_kappa: float = 8.0      # eta
    escalation_sigma: float = 0.35
    p_escalate: float = 0.55
    comp_ref: float = 100.0            # s_ref: component size normaliser
    deg_ref: float = 2.0               # d_ref: repeat-tie degree normaliser
    formation_driver: str = "component"  # "component" | "degree" -- which structural
                                       # quantity sets team-formation search cost.
                                       # "component" (default) makes the observable
                                       # track the slow order parameter; "degree" is
                                       # retained for the sensitivity check.
    output_knee: float = 0.30          # cohesion below which output starts to fall
    output_scale: float = 1.0
    output_sigma: float = 0.10
    output_report_window: int = 10     # program-close outputs are reported as a
                                       # trailing mean over this many cycles (the
                                       # lagging metric's own reporting period)

    # ---- instrument dial 1: consent ----
    consent_c: float = 1.0             # nominal coverage
    consent_arm: str = "mnar"          # "mnar" (centrality-correlated) | "mcar" (uniform)
    consent_beta: float = 1.6          # logistic slope on standardised centrality (MNAR)
    consent_assign_cycle: int = 0      # centrality measured here (end of burn-in),
                                       # then consent is fixed for the recorded run

    # ---- instrument dial 2: proxy error ----
    eps_colour: str = "white"          # "white" | "red" | "systematic"
    eps_severity: str = "mid"          # "low" | "mid" | "high"
    eps_sigma_levels: tuple[float, ...] = (0.4, 0.8, 1.6)   # x observable sd units
    eps_red_phi: float = 0.75          # AR(1) coefficient for red noise
    eps_bias_gain: float = 0.9         # systematic: gain on standardised centrality
    eps_bias_drift: float = 0.004      # systematic: per-cycle instrument drift
    partial_obs_sigma: float = 1.1     # scale of the EXCESS estimation error from
                                       # seeing only part of a team: sigma *
                                       # (1/sqrt(observed dyads) - 1/sqrt(all dyads)),
                                       # which is exactly 0 at full coverage, so c=1.0
                                       # is a clean reference disturbed only by eps

    # ---- instrument dial 3: revocation ----
    revocation: bool = False
    revocation_cycle: int = 140
    revocation_frac: float = 0.40      # fraction of consenters who revoke
    revocation_retroactive: bool = True  # False => static-equivalent comparator

    # ---- instrument dial 4: reward coupling (coverage-and-coupling campaign) ----
    coupling_g: float = 0.0            # 0 | mid | high
    coupling_severed: bool = False     # response computed then discarded (control)
    coupling_gamma: float = 3.0        # proxy-response strength on consented-partner preference
    attrition_kappa: float = 0.02      # Bol-style withdrawal rate per cycle at g=1.
                                       # Withdrawn participants are REPLACED by new
                                       # ones, so the channel raises churn (destroying
                                       # accumulated repeat ties) rather than emptying
                                       # the ecosystem.
    attrition_quantile: float = 0.25   # bottom quantile of measured score withdraws

    # ---- Bifurcation validation: perturb-and-relax probe (relaxation time near the fold) ----
    shock_cycle: int = -1              # recorded cycle at which to shock the tie graph
    shock_frac: float = 0.0            # fraction of repeat ties deleted at the shock

    # ---- bookkeeping ----
    seed: int = 0
    record_series: bool = True

    def a_of_t(self, t: int) -> float:
        if self.scenario == "control":
            return self.a_control
        if t <= self.degrade_start:
            return self.a_high
        span = max(1, self.n_cycles - self.degrade_start)
        frac = (t - self.degrade_start) / span
        return self.a_high + (self.a_low - self.a_high) * frac

    def reengage_prob(self, a: float, z: float) -> float:
        """p_re = A * g(z); g is the cooperative structural feedback term."""
        h2 = self.reengage_hill ** 2
        g = self.reengage_theta + (1.0 - self.reengage_theta) * (z * z) / (h2 + z * z)
        return min(1.0, max(0.0, a * g))

    @property
    def eps_sigma(self) -> float:
        idx = {"low": 0, "mid": 1, "high": 2}[self.eps_severity]
        return self.eps_sigma_levels[idx]

    def to_json(self) -> dict[str, Any]:
        d = asdict(self)
        for k, v in list(d.items()):
            if isinstance(v, tuple):
                d[k] = list(v)
        return d

    @classmethod
    def from_json(cls, d: dict[str, Any]) -> "Config":
        fields = {f for f in cls.__dataclass_fields__}
        kw = {}
        for k, v in d.items():
            if k not in fields:
                continue
            if k in ("team_size_choices", "eps_sigma_levels"):
                v = tuple(v)
            kw[k] = v
        return cls(**kw)


# --------------------------------------------------------------------------
# Small graph helper: windowed cumulative repeat-tie graph
# --------------------------------------------------------------------------
class RepeatTieGraph:
    """Undirected graph of repeat ties, each with an expiry cycle."""

    def __init__(self, n: int, window: int):
        self.n = n
        self.window = window
        self.adj: list[set[int]] = [set() for _ in range(n)]
        self.expiry: dict[tuple[int, int], int] = {}
        self._comp_cache_cycle = -1
        self._comp_of = np.zeros(n, dtype=np.int64)
        self._comp_size = np.zeros(n, dtype=np.int64)   # size of the component of each node

    @staticmethod
    def key(i: int, j: int) -> tuple[int, int]:
        return (i, j) if i < j else (j, i)

    def add(self, i: int, j: int, cycle: int) -> None:
        if i == j:
            return
        k = self.key(i, j)
        self.expiry[k] = cycle + self.window
        self.adj[i].add(j)
        self.adj[j].add(i)
        self._comp_cache_cycle = -1

    def expire(self, cycle: int) -> None:
        dead = [k for k, e in self.expiry.items() if e <= cycle]
        for (i, j) in dead:
            del self.expiry[(i, j)]
            self.adj[i].discard(j)
            self.adj[j].discard(i)
        if dead:
            self._comp_cache_cycle = -1

    @property
    def n_edges(self) -> int:
        return len(self.expiry)

    def components(self, cycle: int) -> tuple[np.ndarray, np.ndarray, list[int]]:
        """Return (comp_of, comp_size_of_node, sizes) with per-cycle caching."""
        if self._comp_cache_cycle == cycle:
            return self._comp_of, self._comp_size, self._sizes
        comp_of = np.full(self.n, -1, dtype=np.int64)
        sizes: list[int] = []
        adj = self.adj
        for start in range(self.n):
            if comp_of[start] != -1:
                continue
            cid = len(sizes)
            stack = [start]
            comp_of[start] = cid
            count = 0
            while stack:
                u = stack.pop()
                count += 1
                for v in adj[u]:
                    if comp_of[v] == -1:
                        comp_of[v] = cid
                        stack.append(v)
            sizes.append(count)
        size_arr = np.asarray(sizes, dtype=np.int64)
        self._comp_of = comp_of
        self._comp_size = size_arr[comp_of]
        self._sizes = sizes
        self._comp_cache_cycle = cycle
        return self._comp_of, self._comp_size, sizes


# --------------------------------------------------------------------------
# The simulation
# --------------------------------------------------------------------------
@dataclass
class RunOutput:
    config: dict[str, Any]
    # per-cycle latent / structural truth
    a_t: list[float] = field(default_factory=list)
    giant_frac: list[float] = field(default_factory=list)
    mean_degree: list[float] = field(default_factory=list)
    susceptibility: list[float] = field(default_factory=list)   # mean finite cluster size
    # per-cycle TRUE observables (full information; not available to the instrument)
    true_formation: list[float] = field(default_factory=list)
    true_escalation: list[float] = field(default_factory=list)
    # per-cycle OBSERVED graph observables: the live repeat-tie graph masked to
    # consented nodes and dyads; nan only when the observed node set is empty
    obs_giant_frac: list[float] = field(default_factory=list)
    obs_mean_degree: list[float] = field(default_factory=list)
    obs_susceptibility: list[float] = field(default_factory=list)
    # per-cycle OBSERVED event observables (consent + proxy error); nan where unobserved
    obs_formation: list[float] = field(default_factory=list)
    obs_escalation: list[float] = field(default_factory=list)
    obs_n_formation: list[int] = field(default_factory=list)
    obs_n_escalation: list[int] = field(default_factory=list)
    # lagging output metric, recorded at program close
    output_close: list[float] = field(default_factory=list)
    output_n: list[int] = field(default_factory=list)
    # coverage accounting (Cencetti input/output split)
    nominal_c: float = 1.0
    achieved_edge_frac: list[float] = field(default_factory=list)
    # Retroactive-revocation second stream: the SAME cycles re-derived from the
    # post-revocation consent set, for every cycle including those published before
    # the revocation.  This is what "revocation propagates through published
    # aggregates" means operationally: the history is re-derived, not appended to.
    obs_formation_rev: list[float] = field(default_factory=list)
    obs_escalation_rev: list[float] = field(default_factory=list)
    # instrument event log
    revocation_applied_cycle: int | None = None
    n_active: list[int] = field(default_factory=list)
    withdrawn_total: int = 0
    # per-program observation records (cycle, team size, observed members,
    # recorded value, true value) -- the Ensign arm's input
    program_log: list[tuple] = field(default_factory=list)
    # raw event tables needed by the auxiliary arms
    events: list[dict[str, Any]] = field(default_factory=list)
    consent_status: list[int] = field(default_factory=list)
    consent_revoked: list[int] = field(default_factory=list)
    observed_degree: list[int] = field(default_factory=list)
    centrality_z: list[float] = field(default_factory=list)

    def series(self) -> dict[str, np.ndarray]:
        return {
            "a_t": np.asarray(self.a_t, float),
            "giant_frac": np.asarray(self.giant_frac, float),
            "mean_degree": np.asarray(self.mean_degree, float),
            "susceptibility": np.asarray(self.susceptibility, float),
            "true_formation": np.asarray(self.true_formation, float),
            "true_escalation": np.asarray(self.true_escalation, float),
            "obs_giant_frac": np.asarray(self.obs_giant_frac, float),
            "obs_mean_degree": np.asarray(self.obs_mean_degree, float),
            "obs_susceptibility": np.asarray(self.obs_susceptibility, float),
            "obs_formation": np.asarray(self.obs_formation, float),
            "obs_escalation": np.asarray(self.obs_escalation, float),
            "obs_formation_rev": np.asarray(self.obs_formation_rev, float),
            "obs_escalation_rev": np.asarray(self.obs_escalation_rev, float),
            "output_close": np.asarray(self.output_close, float),
            "achieved_edge_frac": np.asarray(self.achieved_edge_frac, float),
        }


def _shrink(k_obs: int, m_team: int) -> float:
    """Partial-observation shrinkage of a recorded time-to-team.  Formation latency is
    the span from first to last join; an observer that sees only k of m joins records a
    SHORTER span.  Treating join times as uniform order statistics, the expected
    observed span scales as (k/(k+1)) / (m/(m+1)), which is exactly 1 at full
    observation.  This is the only structural distortion applied to the recorded
    value: the instrument records timestamps, it does not reconstruct the graph."""
    if k_obs >= m_team or m_team <= 1:
        return 1.0
    return (k_obs / (k_obs + 1.0)) / (m_team / (m_team + 1.0))


def _logistic(x: np.ndarray | float) -> np.ndarray | float:
    return 1.0 / (1.0 + np.exp(-x))


def _calibrate_intercept(z: np.ndarray, beta: float, target: float) -> float:
    """Bisect the logistic intercept so mean P(consent) == target."""
    lo, hi = -40.0, 40.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        m = float(np.mean(_logistic(beta * z + mid)))
        if m < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


class Simulation:
    def __init__(self, cfg: Config):
        self.cfg = cfg
        self.rng = make_streams(cfg.seed)
        n = cfg.n_agents
        self.skill = self.rng["intake"].uniform(0.2, 1.0, size=n)
        # friction point 1 (challenge intake / node arrival): every agent is a
        # newcomer until its first participation; `arrived` flips then.
        self.arrived = np.zeros(n, dtype=bool)
        self.experienced = np.zeros(n, dtype=bool)
        self.active = np.ones(n, dtype=bool)          # attrition (coverage-and-coupling campaign) turns these off
        self.collab_count: list[dict[int, int]] = [dict() for _ in range(n)]
        self.last_active = np.full(n, -10**6, dtype=np.int64)
        self.cum_degree = np.zeros(n, dtype=np.int64)  # centrality proxy
        self.measured_score = np.zeros(n, dtype=float)
        self.ties = RepeatTieGraph(n, cfg.tie_window)
        # partner layer: escalation pathways (fixed small set of pathway nodes).
        # Inter-layer edges also live in the temporal window.
        self.n_pathways = 8
        self.pathway_links: list[dict[int, int]] = [dict() for _ in range(self.n_pathways)]
        # consent
        self.consent = np.ones(n, dtype=bool)
        self.consent_final = self.consent
        self.will_revoke = np.asarray([], dtype=np.int64)
        self.revoked = np.zeros(n, dtype=bool)
        self.consent_assigned = False
        # red-noise state
        self._red_state = {"formation": 0.0, "escalation": 0.0}
        # open programs awaiting close (lagging metric)
        self.open_programs: list[dict[str, Any]] = []
        self._out_buffer: list[float] = []
        self._out_cycle_counts: list[int] = []
        self.out = RunOutput(config=cfg.to_json())
        self.out.nominal_c = cfg.consent_c
        # total dyadic team-formation events, and how many were observable
        self._tot_dyads = 0
        self._obs_dyads = 0
        self.n_withdrawn = 0

    # ---------------- consent assignment (dial 1) ----------------
    def _assign_consent(self) -> None:
        cfg = self.cfg
        n = cfg.n_agents
        rng = self.rng["consent"]
        if cfg.consent_c >= 1.0:
            self.consent[:] = True
            self.centrality_z = np.zeros(n)
            self.consent_final = self.consent.copy()
            self.will_revoke = np.asarray([], dtype=np.int64)
            self.consent_assigned = True
            return
        deg = self.cum_degree.astype(float)
        sd = deg.std()
        z = (deg - deg.mean()) / (sd if sd > 1e-9 else 1.0)
        self.centrality_z = z
        if cfg.consent_arm == "mcar":
            p = np.full(n, cfg.consent_c)
        else:  # centrality-correlated => MNAR
            a0 = _calibrate_intercept(z, cfg.consent_beta, cfg.consent_c)
            p = _logistic(cfg.consent_beta * z + a0)
        self.consent = rng.random(n) < p
        self._consent_p = p
        self.consent_assigned = True
        # The revocation set is DRAWN here (at consent assignment) and APPLIED later,
        # at revocation_cycle.  Drawing it up front is what lets the run carry both
        # the as-published stream and its post-revocation re-derivation in one pass,
        # and it keeps the static-equivalent comparator exactly matched: the same
        # agents, the same seed, only the timing differs.
        if cfg.revocation:
            idx = np.flatnonzero(self.consent)
            k = int(round(cfg.revocation_frac * idx.size))
            self.consent_final = self.consent.copy()
            if k > 0 and idx.size:
                z_i = self.centrality_z[idx]
                w = np.asarray(_logistic(0.8 * z_i), float)
                w = w / w.sum()
                self.will_revoke = rng.choice(idx, size=k, replace=False, p=w)
                self.consent_final[self.will_revoke] = False
            else:
                self.will_revoke = np.asarray([], dtype=np.int64)
        else:
            self.consent_final = self.consent
            self.will_revoke = np.asarray([], dtype=np.int64)

    def _apply_revocation(self, t: int) -> None:
        """Apply the pre-drawn revocation set: coverage becomes time-varying MNAR
        (the revoking agents are drawn centrality-weighted from WITHIN the consenting
        population, so who revokes is not who never consented)."""
        if self.will_revoke.size == 0:
            return
        self.consent = self.consent_final.copy()
        self.revoked[self.will_revoke] = True
        self.out.revocation_applied_cycle = t

    # ---------------- proxy error (dial 2) ----------------
    def _proxy_error(self, kind: str, t: int, agent: int, base_sd: float) -> float:
        cfg = self.cfg
        rng = self.rng["proxy"]
        sigma = cfg.eps_sigma * base_sd
        if cfg.eps_colour == "white":
            return float(rng.normal(0.0, sigma))
        if cfg.eps_colour == "red":
            phi = cfg.eps_red_phi
            innov = float(rng.normal(0.0, sigma * math.sqrt(1.0 - phi * phi)))
            self._red_state[kind] = phi * self._red_state[kind] + innov
            return self._red_state[kind]
        # systematic: centrality-linked bias + instrument drift, plus a small white floor
        zc = float(self.centrality_z[agent]) if hasattr(self, "centrality_z") else 0.0
        bias = -cfg.eps_bias_gain * sigma * zc + cfg.eps_bias_drift * sigma * t
        return bias + float(rng.normal(0.0, 0.35 * sigma))

    # ---------------- lifecycle ----------------
    def _shock(self) -> None:
        """Bifurcation validation: delete a random fraction of repeat ties, so the recovery rate of
        the order parameter can be fitted.  The recovery rate going to zero as A
        approaches the fold is the direct signature of a saddle-node bifurcation."""
        rng = self.rng["reengage"]
        keys = list(self.ties.expiry)
        if not keys:
            return
        n_drop = int(round(self.cfg.shock_frac * len(keys)))
        if n_drop <= 0:
            return
        idx = rng.choice(len(keys), size=n_drop, replace=False)
        for m in idx:
            i, j = keys[int(m)]
            del self.ties.expiry[(i, j)]
            self.ties.adj[i].discard(j)
            self.ties.adj[j].discard(i)
        self.ties._comp_cache_cycle = -1

    def _observed_components(self, mask: np.ndarray) -> np.ndarray:
        """Component size, per node, of the subgraph induced on consented nodes by
        consented dyads.  Non-consented nodes get 0 (they are not in the record)."""
        n = self.cfg.n_agents
        adj = self.ties.adj
        seen = np.zeros(n, dtype=bool)
        size_of = np.zeros(n, dtype=np.float64)
        for start in range(n):
            if seen[start] or not mask[start]:
                continue
            stack = [start]
            seen[start] = True
            comp = [start]
            while stack:
                u = stack.pop()
                for v in adj[u]:
                    if mask[v] and not seen[v]:
                        seen[v] = True
                        stack.append(v)
                        comp.append(v)
            size_of[comp] = len(comp)
        return size_of

    def _recycle(self, tabs: int, exclude: list[int]) -> int:
        """Turnover: the longest-inactive agent leaves and a fresh newcomer takes
        its place (same slot index, cleared history).  Keeps the Guimera newcomer
        channel live for the whole run instead of exhausting after ~one window.
        Consent status is a property of the slot and is NOT re-drawn, so nominal
        coverage stays stationary across the run (see RESULTS.md)."""
        cand = np.flatnonzero(self.arrived)
        if cand.size == 0:
            return -1
        if exclude:
            cand = cand[~np.isin(cand, np.asarray(exclude, dtype=np.int64))]
            if cand.size == 0:
                return -1
        i = int(cand[np.argmin(self.last_active[cand])])
        for j in list(self.ties.adj[i]):
            k = self.ties.key(i, j)
            self.ties.expiry.pop(k, None)
            self.ties.adj[j].discard(i)
        self.ties.adj[i].clear()
        self.ties._comp_cache_cycle = -1
        for j in list(self.collab_count[i]):
            self.collab_count[j].pop(i, None)
        self.collab_count[i] = {}
        self.cum_degree[i] = 0
        self.measured_score[i] = 0.0
        self.arrived[i] = False
        self.experienced[i] = False
        self.active[i] = True          # the slot is refilled by a new participant
        return i

    def _pick_team(self, t: int, tabs: int) -> list[int]:
        cfg = self.cfg
        rng = self.rng["formation"]
        m = int(rng.choice(cfg.team_size_choices))
        pool_new = np.flatnonzero((~self.arrived) & self.active)
        pool_exp = np.flatnonzero(self.arrived & self.experienced & self.active)
        pool_any = np.flatnonzero(self.arrived & self.active)
        team: list[int] = []
        for slot in range(m):
            if rng.random() < cfg.p_newcomer:
                cand = -1
                if pool_new.size > 0:
                    cand = int(rng.choice(pool_new))
                    pool_new = pool_new[pool_new != cand]
                elif cfg.turnover:
                    cand = self._recycle(tabs, team)
                if cand >= 0:
                    self.arrived[cand] = True      # friction point 1: node arrival
                    team.append(cand)
                    continue
            # incumbent slot: repeat tie with prob p_repeat_tie (Guimera)
            picked = -1
            if team and rng.random() < cfg.p_repeat_tie:
                anchor = int(rng.choice(np.asarray(team)))
                cands = [j for j in self.ties.adj[anchor]
                         if self.active[j] and j not in team]
                if cfg.coupling_g > 0.0 and cands:
                    picked = self._coupled_pick(cands)
                elif cands:
                    picked = int(rng.choice(np.asarray(cands)))
            if picked < 0:
                base = pool_exp if pool_exp.size > 0 else pool_any
                base = base[~np.isin(base, np.asarray(team, dtype=np.int64))] if team else base
                if base.size == 0:
                    base = pool_any
                if base.size == 0:
                    break
                picked = int(rng.choice(base))
            if picked in team:
                continue
            team.append(picked)
        return team

    def _coupled_pick(self, cands: list[int]) -> int:
        """Dial 4: agents respond to the RECORDED proxy -- they prefer partners whose
        ties are already recorded (consented), because those are the ties the
        published indicator is computed from.  Severed control computes the same
        draw from the same stream, then discards it."""
        cfg = self.cfg
        crng = self.rng["coupling"]
        arr = np.asarray(cands, dtype=np.int64)
        w = 1.0 + cfg.coupling_g * cfg.coupling_gamma * self.consent[arr].astype(float)
        w = w / w.sum()
        coupled_choice = int(crng.choice(arr, p=w))
        uniform_choice = int(crng.choice(arr))
        return uniform_choice if cfg.coupling_severed else coupled_choice

    def _cohesion(self, team: list[int]) -> float:
        if len(team) < 2:
            return 0.0
        pairs = 0
        tied = 0
        for a in range(len(team)):
            for b in range(a + 1, len(team)):
                pairs += 1
                if team[b] in self.ties.adj[team[a]]:
                    tied += 1
        return tied / pairs if pairs else 0.0

    def step(self, t: int) -> None:
        """One program cycle.  `t` is the RECORDED cycle index; it is negative
        during the unrecorded burn-in, during which A is held at a_high and
        nothing is written to the output record."""
        cfg = self.cfg
        record = t >= 0
        tabs = t + cfg.burn_in          # tie-graph / expiry clock
        self.ties.expire(tabs)
        if not self.consent_assigned and t >= cfg.consent_assign_cycle:
            self._assign_consent()
        if (cfg.revocation and cfg.revocation_retroactive
                and t == cfg.revocation_cycle and self.consent_assigned):
            self._apply_revocation(t)
        if t == 0:      # coverage accounting starts with the recorded series
            self._tot_dyads = 0
            self._obs_dyads = 0
        if cfg.shock_frac > 0.0 and t == cfg.shock_cycle:
            self._shock()

        comp_of, comp_size, sizes = self.ties.components(tabs)
        if record:
            A = cfg.a_of_t(t)
        elif cfg.burn_in_a is not None:
            A = cfg.burn_in_a
        else:
            A = cfg.a_control if cfg.scenario == "control" else cfg.a_high
        # structural feedback evaluated at cycle start (the slow mode)
        z_now = 2.0 * self.ties.n_edges / cfg.n_agents
        if (not record) and (t < -cfg.burn_in + cfg.warm_start):
            p_re = min(1.0, A)                 # warm start: feedback bypassed
        else:
            p_re = cfg.reengage_prob(A, z_now)

        adj = self.ties.adj
        deg_full = np.fromiter((len(s) for s in adj), dtype=np.float64,
                               count=cfg.n_agents)
        cons_mask = self.consent
        obs_n_nodes = int(np.count_nonzero(cons_mask))
        if record and obs_n_nodes:
            obs_comp_size = self._observed_components(cons_mask)
            obs_giant_size = float(np.max(obs_comp_size))
            obs_giant_frac = obs_giant_size / obs_n_nodes
            obs_finite = obs_comp_size[(obs_comp_size > 0)
                                       & (obs_comp_size != obs_giant_size)]
            obs_susceptibility = (float(np.mean(obs_finite))
                                  if obs_finite.size else 0.0)
        else:
            obs_giant_frac = math.nan
            obs_susceptibility = math.nan
        # Second observation stream, live only while a drawn-but-unapplied revocation
        # means the as-published mask differs from the post-revocation mask.
        dual = (self.consent_assigned and cfg.revocation
                and not np.array_equal(cons_mask, self.consent_final))
        mask_rev = self.consent_final if dual else cons_mask

        f_true: list[float] = []
        f_obs: list[float] = []
        f_obs_rev: list[float] = []
        e_true: list[float] = []
        e_obs: list[float] = []
        e_obs_rev: list[float] = []

        for _ in range(cfg.programs_per_cycle):
            team = self._pick_team(t, tabs)
            if len(team) < 2:
                continue
            for a in team:
                self.experienced[a] = True
                self.last_active[a] = tabs

            # --- friction point 2: team formation -> formation latency -------
            # Search cost falls with the size of the repeat-tie component the team
            # can recruit from; averaged over the team so the observable tracks the
            # slow structural mode rather than one anchor's idiosyncrasy.
            shock = float(self.rng["formation"].normal(0.0, cfg.formation_sigma))
            tarr = np.asarray(team)
            if cfg.formation_driver == "component":
                avail = float(np.mean(comp_size[tarr])) / cfg.comp_ref
            else:
                avail = float(np.mean(deg_full[tarr])) / cfg.deg_ref
            lat = max(0.1, cfg.formation_base + cfg.formation_kappa / (1.0 + avail) + shock)
            f_true.append(lat)

            # Observation: the instrument sees the formation event only through
            # consented dyads, and reconstructs the latency from the consented
            # members' structure alone -- so under MNAR the recorded value is
            # structurally biased, not merely noisier.
            cons = [a for a in team if self.consent[a]]
            n_obs_dyads = len(cons) * (len(cons) - 1) // 2
            n_dyads = len(team) * (len(team) - 1) // 2
            self._tot_dyads += n_dyads
            self._obs_dyads += n_obs_dyads
            m_team = len(team)
            if n_obs_dyads >= 1:
                excess = (1.0 / math.sqrt(n_obs_dyads)) - (1.0 / math.sqrt(n_dyads))
                est_noise = float(self.rng["proxy"].normal(
                    0.0, cfg.partial_obs_sigma * max(0.0, excess)))
                err = self._proxy_error("formation", max(0, t), cons[0],
                                        cfg.formation_kappa * 0.12)
                v = lat * _shrink(len(cons), m_team) + est_noise + err
                f_obs.append(v)
                if record:
                    self.out.program_log.append((t, m_team, len(cons), v, lat))
                if dual:
                    k_r = sum(1 for a in team if mask_rev[a])
                    if k_r >= 2:
                        f_obs_rev.append(lat * _shrink(k_r, m_team)
                                         + est_noise + err)
            elif dual:
                k_r = sum(1 for a in team if mask_rev[a])
                if k_r >= 2:
                    f_obs_rev.append(lat * _shrink(k_r, m_team))

            # --- friction point 3: sprint progression -----------------------
            coh = self._cohesion(team)
            prng = self.rng["progress"]
            stall = float(prng.random() < (0.55 - 0.45 * coh))
            iterations = 1 + int(prng.poisson(1.2 + 2.0 * coh))

            # --- friction point 4: evaluation (edge weighting) ---------------
            concord = float(prng.normal(0.6 + 0.3 * coh, 0.12))

            # --- friction point 5: escalation (inter-layer edge) -------------
            erng = self.rng["escalate"]
            if erng.random() < cfg.p_escalate:
                eshock = float(erng.normal(0.0, cfg.escalation_sigma))
                reach = self._pathway_reach(team, comp_of, tabs)
                elat = max(0.1, cfg.escalation_base
                           + cfg.escalation_kappa / (1.0 + reach) + eshock)
                e_true.append(elat)
                esc_agent = team[int(erng.integers(len(team)))]
                if self.consent[esc_agent]:
                    err = self._proxy_error("escalation", max(0, t), esc_agent,
                                            cfg.escalation_kappa * 0.12)
                    e_obs.append(elat + err)
                    if dual and mask_rev[esc_agent]:
                        e_obs_rev.append(elat + err)
                elif dual and mask_rev[esc_agent]:
                    e_obs_rev.append(elat)
                pw = int(erng.integers(self.n_pathways))
                self.pathway_links[pw][esc_agent] = tabs + cfg.tie_window

            # --- friction point 6: partner re-engagement -> repeat ties ------
            rrng = self.rng["reengage"]
            for a in range(len(team)):
                for b in range(a + 1, len(team)):
                    i, j = team[a], team[b]
                    self.collab_count[i][j] = self.collab_count[i].get(j, 0) + 1
                    self.collab_count[j][i] = self.collab_count[j].get(i, 0) + 1
                    self.cum_degree[i] += 1
                    self.cum_degree[j] += 1
                    if rrng.random() < p_re:
                        self.ties.add(i, j, tabs)

            # program stays open; the lagging output metric lands at close
            self.open_programs.append({
                "close": tabs + cfg.program_duration,
                "team": team, "cohesion": coh, "stall": stall,
                "iterations": iterations, "concord": concord,
            })
            for a in team:
                self.measured_score[a] = 0.85 * self.measured_score[a] + 0.15 * (
                    iterations / 4.0 + coh)

        # --- lagging output metric: programs closing this cycle --------------
        outs: list[float] = []
        still: list[dict[str, Any]] = []
        for pr in self.open_programs:
            if pr["close"] <= tabs:
                coh = pr["cohesion"]
                # saturating sensitivity: output holds up until cohesion drops
                # below the knee, then falls -- programs absorb coordination loss
                # before it shows in output (this is why a lead time can exist)
                q = min(1.0, coh / self.cfg.output_knee)
                val = (self.cfg.output_scale * (0.35 + 0.65 * q)
                       * (1.0 - 0.25 * pr["stall"])
                       + float(self.rng["progress"].normal(0.0, self.cfg.output_sigma)))
                outs.append(val)
            else:
                still.append(pr)
        self.open_programs = still

        # --- coverage-and-coupling campaign attrition channel (Bol et al. 2018) ------------------------
        if cfg.coupling_g > 0.0:
            self._attrition(t)

        # --- record -----------------------------------------------------------
        if not record:
            return
        gc = max(sizes) if sizes else 0
        finite = [s for s in sizes if s != gc] or [0]
        self.out.a_t.append(A)
        self.out.giant_frac.append(gc / cfg.n_agents)
        self.out.mean_degree.append(2.0 * self.ties.n_edges / cfg.n_agents)
        obs_n_edges = sum(1 for i, j in self.ties.expiry
                          if cons_mask[i] and cons_mask[j])
        self.out.obs_giant_frac.append(obs_giant_frac)
        self.out.obs_mean_degree.append(2.0 * obs_n_edges / obs_n_nodes
                                        if obs_n_nodes else math.nan)
        self.out.obs_susceptibility.append(obs_susceptibility)
        self.out.susceptibility.append(float(np.mean(np.asarray(finite, float) ** 2)
                                             / max(1.0, np.mean(finite))))
        self.out.true_formation.append(float(np.mean(f_true)) if f_true else math.nan)
        self.out.true_escalation.append(float(np.mean(e_true)) if e_true else math.nan)
        self.out.obs_formation.append(float(np.mean(f_obs)) if f_obs else math.nan)
        self.out.obs_escalation.append(float(np.mean(e_obs)) if e_obs else math.nan)
        self.out.obs_formation_rev.append(float(np.mean(f_obs_rev))
                                          if f_obs_rev else
                                          (float(np.mean(f_obs)) if f_obs else math.nan))
        self.out.obs_escalation_rev.append(float(np.mean(e_obs_rev))
                                           if e_obs_rev else
                                           (float(np.mean(e_obs)) if e_obs else math.nan))
        self.out.obs_n_formation.append(len(f_obs))
        self.out.obs_n_escalation.append(len(e_obs))
        # The lagging metric is what a program-output report actually publishes:
        # the mean over the reporting period, not one cycle's programs.
        self._out_buffer.extend(outs)
        self._out_cycle_counts.append(len(outs))
        if len(self._out_cycle_counts) > cfg.output_report_window:
            drop = self._out_cycle_counts.pop(0)
            if drop:
                del self._out_buffer[:drop]
        self.out.output_close.append(float(np.mean(self._out_buffer))
                                     if self._out_buffer else math.nan)
        self.out.output_n.append(len(self._out_buffer))
        self.out.achieved_edge_frac.append(
            self._obs_dyads / self._tot_dyads if self._tot_dyads else math.nan)
        self.out.n_active.append(int(self.active.sum()))

    def _pathway_reach(self, team: list[int], comp_of: np.ndarray, t: int) -> float:
        """Number of escalation pathways reachable from the team through the
        current repeat-tie structure -- the inter-layer routing capacity."""
        comps = {int(comp_of[a]) for a in team}
        reach = 0
        for pw in range(self.n_pathways):
            links = self.pathway_links[pw]
            for member, exp in list(links.items()):
                if exp <= t:
                    del links[member]
                    continue
                if int(comp_of[member]) in comps:
                    reach += 1
                    break
        return float(reach)

    def _attrition(self, t: int) -> None:
        """Bol-style withdrawal: agents measured in the bottom quantile of the
        recorded score withdraw with probability g * kappa."""
        cfg = self.cfg
        crng = self.rng["coupling"]
        act = np.flatnonzero(self.active & self.experienced)
        if act.size < 10:
            return
        thr = float(np.quantile(self.measured_score[act], cfg.attrition_quantile))
        low = act[self.measured_score[act] <= thr]
        if low.size == 0:
            return
        draws = crng.random(low.size)
        leave = low[draws < cfg.coupling_g * cfg.attrition_kappa]
        if cfg.coupling_severed:
            return                    # stream consumed, effect discarded
        self.active[leave] = False
        self.n_withdrawn += int(leave.size)

    def run(self) -> RunOutput:
        for t in range(-self.cfg.burn_in, self.cfg.n_cycles):
            self.step(t)
        self.out.withdrawn_total = int(self.n_withdrawn)
        self.out.consent_status = self.consent.astype(int).tolist()
        self.out.consent_revoked = self.revoked.astype(int).tolist()
        # observed degree = degree in the observed (consented-dyad) repeat-tie graph
        obsdeg = np.zeros(self.cfg.n_agents, dtype=np.int64)
        for (i, j) in self.ties.expiry:
            if self.consent[i] and self.consent[j]:
                obsdeg[i] += 1
                obsdeg[j] += 1
        self.out.observed_degree = obsdeg.tolist()
        self.out.centrality_z = (self.centrality_z.tolist()
                                 if hasattr(self, "centrality_z")
                                 else [0.0] * self.cfg.n_agents)
        return self.out


def run_config(cfg: Config) -> RunOutput:
    return Simulation(cfg).run()


def run_from_dict(d: dict[str, Any]) -> RunOutput:
    return run_config(Config.from_json(d))


if __name__ == "__main__":
    import time
    for scen in ("degrade", "control"):
        c = Config(scenario=scen, seed=1)
        t0 = time.time()
        o = run_config(c)
        s = o.series()
        print(f"{scen}: {time.time()-t0:.2f}s  giant[40]={s['giant_frac'][40]:.3f} "
              f"giant[-1]={s['giant_frac'][-1]:.3f}  z[40]={s['mean_degree'][40]:.2f} "
              f"z[-1]={s['mean_degree'][-1]:.2f} "
              f"form[40]={s['true_formation'][40]:.2f} form[-1]={s['true_formation'][-1]:.2f} "
              f"out[50]={s['output_close'][50]:.3f} out[-1]={s['output_close'][-1]:.3f}")
    print(json.dumps({"ok": True}))
