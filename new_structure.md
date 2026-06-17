## Central Research Question

To what extent do machine learning-based and classical kappa weighting estimators
of the LATE implicitly target different subpopulations, and what do their outcome
weights reveal about covariate balance and estimator reliability across empirical
applications?


### Sub-questions I want to answer
**RQ1.** What are the structural properties of outcome weights for kappa-based LATE (ESS- here there two formulas that defined that, negative weight share, sum-to-zero as translation invariance
criterion), and how do they differ from those of DML-based estimators independently
of any specific dataset?
Importnat to bridge teh gap between what he defines as outcome wieghts statistics and what I define in my function


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
1. Derive closed-form outcome weights for τ̂ᵤ in the Knaus PIVE framewor.
2. Show translation invariance for all estimators. This should be done using the rerun check which means fitting the etsimator twice using the same random seed botj times and then check the fifference is the same as log(100) times the weight
3. Clarify the distinction between Abadie's kappa weights (identification objects)
   and outcome weights in the PIVE sense (ωᵢ such that τ̂ = ΣωᵢYᵢ).
4. Apply Love maybe the package offers something unique for the binary variables he meant when it comes to depicting plots and ESS diagnostics to kappa estimators for the first time,
   using the same pipeline as Knaus (2024).
5. Compare kappa estimators (τ̂ᵤᵐˡ, τ̂ᵤᶜᵇ, τ̂ₐ,₁₀) with DML Wald-AIPW across
   three empirical applications, using multiple ML learners for the nuisance
   parameters.
6. Apply the DML estimators using Knaus outcome weights package. And make a full diagnostics section meaning really look sometimes at the smoother matrix and maybe look also at the descritptives or so.  
7. Implement the Wald AIPW estimation with double machine learning using XGBOOST and linear regression and the ranger function of the package and try there to set up a tuning section which tunes the parameters for the xgboost in order to get clean results for the DMl estimator section. Moreover show what the difference is between what he implemented there in his developmemt package and what is done with the outcome weights package
8. Check descriptives for all of the three framewoks, things such as, how good the first stage, is how propensity score differs firn really the instruement values, whetehr teh isntrument and everything is binary and there are no missing or so, make it clean so that no errros appear
9. New task from Knaus whihc he summarised as 2. Eine sehr spannende Erweiterung wäre noch Wald-AIPW und PLR-IV mit logit pscore und OLS outcome regression zu implementieren, quasi als Zwischenstufe zwischen den parametrischen Methoden im Wooldridge Paper und den DML Methoden. Das ist nicht im OutcomeWeights Paket abgedeckt und wäre gerade deswegen ideal, um auch ein bisschen zu Coden. Am einfachsten wäre es vermutlich, das über das DoubleML package zu machen. Ich hänge dir mal an, wie man es ohne Instrument innerhalb von DoubeML implementieren würde. Ich denke es würde dann darauf hinauslaufen diese Funktion https://github.com/MCKnaus/OutcomeWeights/blob/0c94f940b04c14d0247b46842af37752e306b79e/R/DoubleML.R#L191 kompatibel zu machen mit lrn("regr.lm"). Wenn du das sauber hinbekommst (ich würde dir bei Problemen auch helfen) hättest du deutlich mehr erreicht, als wenn du die instrumental_forest Funktionen einfach anwendest und es würde sich exzellent in die übergeordnete Fragestellung einfügen. Was meinst du?


---


### Chapter 2 — Econometric Framework (8–10 pages) 🔲

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

**Section 2.3 — Kappa-based LATE estimators**
- The five estimators: τ̂ᵤ (Uysal 2011), τ̂ₐ,₁₀ (Abadie & Cattaneo 2018),
  unnormalized τ̂ₐ, τ̂ₜ (= τ̂ₐ,₁, Frölich/Tan), τ̂ₐ,₀
- Normalized vs. unnormalized: what the distinction means mechanically
- Propensity score estimation: MLE logit (τ̂ᵤᵐˡ) vs. CBPS (τ̂ᵤᶜᵇ);
  Proposition 3.5: with CBPS all normalized estimators coincide

**Section 2.4 — Why normalization matters**
- Definition TI (translation invariance): τ̂(Y, W) = τ̂(Y+k, W) for all k
- Proposition 3.2 (SUW 2025): τ̂ᵤ and τ̂ₐ,₁₀ pass; τ̂ₐ, τ̂ₜ, τ̂ₐ,₀ fail
- Definition SE (scale equivariance): brief statement, linked to log-unit sensitivity
- Concrete example: cents vs. dollars failure from Table 2 of SUW 202, here especially later mentioning that there are 2 methods how to measure it but that I am keen on method A as this is even proposed by the professor & this is the way it is gonna be implemented in the code

**Section 2.5 — DML and Wald-AIPW**
- DML framework (Chernozhukov et al. 2018): PLR-IV and Wald-AIPW estimator
- -Maybe set it a bit apart like how PLR-IV does not always measure LATE
- The difference between LATE and the constant structural treatment theta and when
  they are the same
- Two nuisance parameters: E[Y|Z, X] and E[D|Z, X], estimated via K-fold
  cross-fitting
- - For XGBoost, hyperparameters are selected by inner cross-validation.
- - The tuning criterion is predictive nuisance loss, not the causal estimand.


**Section 2.6 — PIVE framework and outcome weights (Knaus 2024)**
- Any estimator fitting the pseudo-IV structure: τ̂ = Σᵢ ωᵢYᵢ
- The two-step: (i) form pseudo-instrument Z̃ and transformation matrix T;
  (ii) ω' = (Z̃'D̃)⁻¹Z̃'T
- "Fully normalized" in Knaus (2024): Σ_{D=1} ωᵢ = +1, Σ_{D=0} ωᵢ = −1
  (Table 5)
**Section 2.7 — Covariate balance diagnostics**
- Standardized Mean Difference (SMD): |X̄ₜᵣₑₐₜₑ_ₖ − X̄_cₒₙₜᵣₒₗ_ₖ| / SD(Xₖ),
  computed with outcome weights ωᵢ
- Love plots: one dot per covariate, unadjusted vs. weighted SMD; threshold at
  |SMD| ≤ 0.1
- Effective Sample Size (ESS): ESS = 1 / Σᵢ ωᵢ², also to take into account are the modified ESS mentioned in Zubizaretta
- Negative weight share: % of observations with ωᵢ < 0

---


### Chapter 3 — Connecting the Frameworks 


**Section 3.1 — Kappa weights vs. outcome weights** ✅
- Distinction between κᵢ (identification objects, Abadie 2003) and ωᵢ
  (sample-level PIVE weights, Knaus 2024) clearly established
- Subsection 3.1.1 (Two notions of normalization): SUW normalization (estimator
  construction) vs. Knaus normalization
- With regard to why do we need outcome weights. Love plots use ωᵢ, not κᵢ — stated explicitly

**Section 3.2 — PIVE framework** 
- Proposition 1 (Knaus 2024) stated formally with label 
- Diagonal T matrix structure for kappa estimators established — no smoother
  condition needed; existence of ωᵢ is purely algebraic


**Section 3.3 — Outcome weights for τ̂ᵤ** 
- Full derivation: normalized IPW contrast → diagonal T^u → closed-form ωᵢᵘ
  (eq:u_omega_scalar, boxed)
- Three normalization conditions verified algebraically:
  Σωᵢ = 0 (via equal-mass property), ΣωᵢDᵢ = 1, Σωᵢ(1−Dᵢ) = −1
- Remark (rem:hajek_contrast): why τ̂ₐ,₁ fails where τ̂ᵤ succeeds, the Hájek
  normalization is the single algebraic step that determines translation
  invariance; connects to Appendix E derivation


**Section 3.4 — Normalization classification and summary** 
- Subsection 3.4.1: how SUW normalization and Knaus normalization align for the
  kappa family — co-occurrence is algebraic, not a general theorem
- Subsection 3.4.2: Table 1 (tab:normalization) — all six estimators classified
  by SUW norm., Σωᵢ=0, ΣωᵢDᵢ, Σωᵢ(1−Dᵢ), Knaus class


**Section 3.5 — From Derivation to Diagnostics: Computational Implementation** 
- kappa_outcome_weights(Z, D, p) described: returns all five ωᵢ vectors in
  closed form, no numerical optimisation
- check_weight_identity() and weight_diag() described as companion functions
- Pipeline: propensity score → weights → verify identity → Love plots 
---


### Chapter 4 — Empirical Application: Angrist (1990) Vietnam Draft Lottery

1. Data and the design

- The draft Lottery -> the instrument

- Smaple, variables and design Diagnostics

2. Point Estimates and Replication

3. Double Machine Learning Comparison

- Comparing the double machine leanring estimators of his outcome weights package with the ones in the outcome package development function

- Maybe also implementing the last task extension he gave me to get the smoother matrix somehow as well using the linear or the logistic regression.

4.Outcome Weights Diagnostics and Covariate Balance

- Here trying to make a mix of my outcome weight diagnostic and his summary that he gave me using just summary of his package

- Make translation invariance visible and see if values change why this is the case 

5. Love plots

6. Conclusion and on Angrist Vietnam Framewrok







### Chapter 5 — Empirical Applications: Card (1995) and Angrist & Evans (1998)

1. Data 

- College Proximity - the instrument

- Treatment Definintions and outcomes

- Covariate Specifications

2. Point estimates and translation invariance




### Chapter 6 - Angrist and Evans Labour Supply

1. Data

- For now find the sample size: It has 290 000 obs which is way to much and would cause the smoother matrix to crash, thinking off how can I thoughtfully choose from the sample maybe a small sample, 3-5k lets see what happens and who do I choose.


### Chapter 7 — Discussion 

1. Creating a cross Application summary, especailly being key on a comaprative design diagnostic table across all three applications.

2. What does the outcome weights lens add

3. DML learner comparison
- especially what my new implementation gave me based on what knaus said, to be figured out.


### Chapter 8 — Conclusion (1–2 pages)