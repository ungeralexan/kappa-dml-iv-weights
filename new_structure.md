## Central Research Question

To what extent do machine learning-based and classical kappa weighting estimators
of the LATE implicitly target different subpopulations, and what do their outcome
weights reveal about covariate balance and estimator reliability across empirical
applications?


### Sub-questions I want to answer
**RQ1.** What are the structural properties of outcome weights for kappa-based
LATE estimators, and how do they differ from those of DML-based estimators
independently of any specific dataset?

**RQ2.** When applied to empirical datasets, how do outcome weights of DML-based
estimators (Wald-AIPW with cross-fitted random forests) and kappa estimators compare
in terms of covariate balance (Love plots / SMDs) and negative weights patterns?

**RQ3.** Can outcome-weight diagnostics guide practitioners toward more robust
estimator choices, and does the outcome-weights lens explain divergences between
classical kappa, normalized kappa, and DML-based IV estimators?

## Thesis structure (Gliederung)

### Chapter 1 — Introduction 
**Section 1.1 — Motivation**
Hook: the same dataset yields wildly different LATE estimates depending on which
estimator is used — not because of different assumptions, but because of how the
outcome variable is coded. Unnormalized kappa estimators violate translation
invariance: adding a constant to every outcome changes the treatment effect estimate.
This is the practical failure mode that motivates the thesis.

**Section 1.2 — Background**

**Section 1.3 — Research gap**
Knaus (2024) introduces the PIVE framework and derives outcome weights ωᵢ such that
τ̂ = Σᵢ ωᵢYᵢ for DML/GRF estimators, enabling covariate balance diagnostics via
Love plots. Appendix A.4 of Knaus (2024) sketches the same derivation for kappa
estimators but does not apply it empirically. This thesis fills that gap.


**Section 1.4 — Contribution**
1. Derive closed-form outcome weights for τ̂ᵤ in the Knaus PIVE framework.
2. Separate three related concepts that are easy to confuse: SUW estimator
   normalization, Knaus outcome-weight normalization, and translation-invariance
   rerun checks.
3. Clarify the distinction between Abadie's kappa weights (identification objects)
   and outcome weights in the PIVE sense (ωᵢ such that τ̂ = ΣωᵢYᵢ).
4. Apply outcome-weight diagnostics, including Love plots, SMDs, negative weights,
   and ESS measures, to kappa estimators using the same diagnostic logic as
   Knaus (2024).
5. Compare kappa estimators (τ̂ᵤᵐˡ, τ̂ᵤᶜᵇ, τ̂ₐ,₁₀) with DML Wald-AIPW across
   three empirical applications, using multiple ML learners for the nuisance
   parameters.
6. Implement and compare DML learner variants for Wald-AIPW and PLR-IV,
   including linear/logit baselines, Ranger, and tuned/untuned XGBoost.
7. Use Method-A rerun checks to test translation invariance of the implemented
   estimators, controlling seeds, folds, learners, and tuning rules.
8. Provide clean descriptive and design diagnostics for each empirical application
   before interpreting outcome-weight results.


---


### Chapter 2 The Econometric Framework 

**Section 2.1 — IV, LATE, and compliers**
- Potential outcomes notation: Yᵢ(0), Yᵢ(1), Dᵢ(0), Dᵢ(1)
- Four compliance types (AIR 1996): always-takers, never-takers, compliers, defiers
- IV assumptions (i)–(iv): conditional independence, exclusion restriction,
  first stage / overlap, monotonicity (I still have to cites those)
- LATE definition: τᴸᴬᵀᴱ = E[Y₁ − Y₀ | D₁ > D₀]
- Properties of a poor instrument?
- Why 2SLS may not recover LATE under heterogeneous effects (one-sentence reference
  to Blandhol et al. 2022)


**Section 2.2 — Abadie's kappa theorem**
- Lemma 2.1 (Abadie 2003, restated in SUW 2025 notation): the three weights
  κ, κ₁, κ₀ and their cell-by-cell values (Table 1 of SUW 2025)
- Parts (a), (b), (c) of the kappa theorem: any complier moment is identified
- Remark 2.2: E(κ) = E(κ₁) = E(κ₀) = P(D₁ > D₀) in population; why they
  diverge in finite samples


### Chapter 3 The Estimators


**Section 3.1 — Kappa-based LATE estimators**
- The five estimators: τ̂ᵤ (Uysal 2011), τ̂ₐ,₁₀ (Abadie & Cattaneo 2018),
  unnormalized τ̂ₐ, τ̂ₜ (= τ̂ₐ,₁, Frölich/Tan), τ̂ₐ,₀
- Normalized vs. unnormalized: what the distinction means mechanically
- Propensity score estimation: MLE logit (τ̂ᵤᵐˡ) vs. CBPS (τ̂ᵤᶜᵇ);
  Proposition 3.5: with CBPS all normalized estimators coincide

**Section 3.2 — DML-based IV estimators**
- DML framework (Chernozhukov et al. 2018): PLR-IV and Wald-AIPW estimator
- Set PLR-IV apart clearly: it is a partially linear IV estimator and does not
  automatically target the same LATE object unless the structural assumptions
  justify that interpretation
- The difference between LATE and the constant structural treatment theta and when
  they are the same
- Two nuisance parameters: E[Y|Z, X] and E[D|Z, X], estimated via K-fold
  cross-fitting
- For XGBoost, hyperparameters are selected by inner cross-validation
- The tuning criterion is predictive nuisance loss, not the causal estimand


### Chapter 4  Connecting the Frameworks
**Section 4.1 — Outcome weights and Knaus normalization**
- Introduce outcome weights conceptually: τ̂ = Σᵢ ωᵢYᵢ
- The distinction between Abadie’s Ki and and Knaus’s ωᵢ
  - ki identification weights derived from the population complier representation;
  - ωᵢ final sample-level weights that reproduce the fitted estimator numerically.
- The two-step: (i) form pseudo-instrument Z̃ and transformation matrix T;
  (ii) ω' = (Z̃'D̃)⁻¹Z̃'T
- Knaus classes:
  - translation/scale-normalized: Σᵢωᵢ = 0
  - fully normalized: Σ_{D=1}ωᵢ = +1 and Σ_{D=0}ωᵢ = −1
- Full derivation: normalized IPW contrast → diagonal T^u → closed-form ωᵢᵘ
  (eq:u_omega_scalar, boxed)
- Three normalization conditions verified algebraically:
  Σωᵢ = 0 (via equal-mass property), ΣωᵢDᵢ = 1, Σωᵢ(1−Dᵢ) = −1
- Remark (rem:hajek_contrast): why τ̂ₐ,₁ fails where τ̂ᵤ succeeds, the Hájek
  normalization is the single algebraic step that determines translation
  invariance if this could go into Appendix
- **From Derivation to Diagnostics: Computational Implementation** 
- kappa_outcome_weights(Z, D, p) described: returns all five ωᵢ vectors in
  closed form, no numerical optimisation
- check_weight_identity() and weight_diag() described as companion functions
- Pipeline: propensity score → weights → verify identity → Love plots 



**Section 4.2 — Translation and scale invariance**
- Definition TI: τ̂(Y, W) = τ̂(Y+k, W) for all k (translation invariance) (Also think about the binary recoding case in the last framework)
- Outcome-weight proof: τ̂(Y+k) − τ̂(Y) = kΣᵢωᵢ
- Proposition 3.2 (SUW 2025): τ̂ᵤ and τ̂ₐ,₁₀ pass; τ̂ₐ, τ̂ₜ, τ̂ₐ,₀ fail
- Explain scale equivariance differently Definition SE / scale issue: brief statement, linked to log-unit sensitivity
- Concrete example: cents vs. dollars is an additive shift after logs
- Separate clearly:
  - Method B = frozen-weight algebraic check using extracted ωᵢ
  - Method A = full rerun of the implemented estimator on Y+k
- Emphasize that Method A is the later empirical contribution because it includes
  nuisance training, cross-fitting, tuning, and algorithmic randomness

**Section 4.3 — Outcome weights diagnostics**
- Covariate Balance: Standardized Mean Difference (SMD): |X̄ₜᵣₑₐₜₑ_ₖ − X̄_cₒₙₜᵣₒₗ_ₖ| / SD(Xₖ),
  computed with outcome weights ωᵢ
- Love plots: one dot per covariate, unadjusted vs. weighted SMD; threshold at
  |SMD| ≤ 0.1
- Effective Sample Size (ESS): standard Kish-style ESS and the modified ESS used
  in the thesis diagnostics
- Negative weight share: % of observations with ωᵢ < 0
- Extreme weights concentration
 - maximum absolute weight;
 - upper quantiles of wi
 - top -1% share of absolute weight mass 

---


### Chapter 5 — Empirical Application: Angrist (1990) Vietnam Draft Lottery

1. Data and the design

- Draft lottery instrument
- Sample, variables, and design diagnostics

2. Point Estimates and Replication

3. Double Machine Learning Comparison

- Compare DML estimators from the OutcomeWeights implementation with the
  manually implemented learner variants where relevant
- Include the intermediate linear/logit implementation if it strengthens the
  bridge between kappa estimators and flexible DML estimators

4. Outcome Weights Diagnostics and Covariate Balance

- Combine Knaus-style summary statistics with the thesis-specific diagnostics
- Make translation-invariance sensitivity visible where the outcome coding changes

5. Love plots

6. Short application conclusion





### Chapter 6 — Empirical Application: Card (1995)

1. Data 

- College Proximity - the instrument

- Treatment definitions and outcomes

- Covariate Specifications

2. Point estimates and translation invariance




### Chapter 7 — Empirical Application: Angrist & Evans (1998) Childbearing

1. Data and sample construction

- Instrument: same-sex composition of the first two children
- Treatment: more than two children
- Outcomes:
  - labor supply / worked last year (`workedm`)
  - log income (`lincomem`)
- Two relevant analysis samples:
  - labor/binary-outcome sample
  - positive-income subsample for log income
- Diagnostic subsamples: honest stratified 3,000-observation draws, one for labor
  and one for income, designed to preserve the IV structure without seed searching

2. Descriptive diagnostics

3. Kappa replication and SUW recodings

- Replicate the normalized and unnormalized kappa estimators from SUW
- Binary outcome recodings:
  - `workedm`
  - additive translation check where relevant
  - display convention for “did not work” kept separate from the actual analysis
- Income recodings:
  - log income in different monetary units


4. DML and learner comparison

- DML smoother / GRF-style headline estimators: PLR-IV and Wald-AIPW
- Learner comparison:
  - linear/logit baseline
  - Ranger
  - XGBoost imported from separate tuning scripts
- PLR-IV and Wald-AIPW shown separately to avoid mixing estimator families
- XGBoost tuning:
  - untuned baseline kept fixed
  - tuned version selected by inner CV on nuisance-prediction loss
  - weights extracted where available and flagged when identity checks fail


5. Translation-invariance rerun check

- Method B: frozen-weight algebraic prediction from extracted outcome weights
- Method A: full rerun on shifted outcomes with same seed, folds, learner, and tuning rule
- Labor: additive shift of `workedm`
- Income: log-unit shift such as `log(100)`
- XGBoost Method A:
  - original tuned estimates imported from saved tuning exports
  - shifted outcome rerun with the same nested tuning procedure
- Interpretation of results belongs in the findings chapter; this chapter reports
  the diagnostic structure and main checks

### Chapter 8 — Discussion 

1. Cross-application summary, with a comparative design-diagnostic table across
   all empirical applications

2. What does the outcome weights lens add

3. DML learner comparison
- Discuss what the additional linear/Ranger/XGBoost implementations add relative
  to the package-based diagnostics


### Chapter 8 — Conclusion (1–2 pages)
