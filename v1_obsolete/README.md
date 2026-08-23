# Ring Vaccination Operating Characteristics Model

## Methodology and Simulation Framework

## 1. Overview

This model evaluates the operating characteristics of pragmatic ring
vaccination study designs under uncertainty in:

-   vaccine biological effectiveness;
-   transmission intensity;
-   between-ring heterogeneity in transmission risk;
-   number of vaccinated rings;
-   allocation ratio between immediate and delayed vaccination;
-   vaccination delay.

The model is designed to quantify the trade-off between:

1.  **statistical efficiency**: ability to estimate vaccine
    effectiveness (VE) precisely; and
2.  **ethical and operational efficiency**: maximising cases prevented
    while minimising unnecessary delay to vaccination.

Rather than evaluating only conventional trial power, the framework
explicitly considers the practical decision problem faced during
outbreak response: how much randomisation or delay is justified given
the competing goals of learning and preventing disease.

------------------------------------------------------------------------

# 2. Conceptual model

The simulation separates the epidemic process from the intervention
design.

This allows multiple vaccination strategies to be evaluated against the
same underlying epidemic realisations.

The structure is:

    Generate epidemic scenarios
            |
            v
    Assign latent ring-level transmission risk
            |
            v
    Create epidemic cache
            |
            v
    Apply vaccination design
            |
            +--> allocation ratio
            +--> vaccination delay
            +--> vaccine effectiveness
            |
            v
    Estimate VE and public health outcomes

This separation reduces Monte Carlo noise when comparing alternative
designs.

------------------------------------------------------------------------

# 3. Epidemic simulation layer

## 3.1 Ring structures

Each simulation consists of a defined number of vaccination rings.

The number of rings is varied as a design parameter.

Current operating range:

-   400 rings
-   600 rings
-   1000 rings
-   1500 rings

Each ring contains:

-   index case;
-   contacts;
-   contacts of contacts (CoCs), where included;
-   exposure structure;
-   baseline transmission risk.

------------------------------------------------------------------------

## 3.2 Transmission probability

The model represents transmission opportunity using the secondary attack
rate (SAR).

The primary SAR scenarios are:

  Parameter               Values
  ----------- ------------------
  SAR           0.02, 0.05, 0.10

These represent low, moderate, and higher transmission settings.

The mean SAR is preserved when introducing heterogeneity.

------------------------------------------------------------------------

# 4. Between-ring SAR heterogeneity

A key extension of the model is inclusion of latent variation in SAR
between rings.

A single average SAR may hide substantial differences between outbreak
clusters. Some rings may have no secondary infections, while others may
have several.

The model therefore includes a coefficient of variation (CV) parameter:

  Parameter           Values
  ----------- --------------
  SAR CV        0, 0.25, 0.5

Conceptually:

-   SAR CV = 0 assumes all rings have identical transmission potential.
-   SAR CV \> 0 allows some rings to be intrinsically higher-risk than
    others.

The purpose is not to change expected transmission, but to assess how
hidden heterogeneity affects trial inference.

Expected effect:

    Higher SAR CV
          |
          v
    greater outcome variability
          |
          v
    larger uncertainty in VE estimate

------------------------------------------------------------------------

# 5. Vaccination design layer

After epidemic structures are generated, vaccination strategies are
applied.

## 5.1 Allocation ratio

The model allows different proportions of rings to receive immediate
versus delayed vaccination.

Examples:

-   1:1
-   1.5:1
-   2:1
-   3:1

The allocation ratio affects:

-   ethical trade-off;
-   number of cases prevented;
-   statistical precision.

------------------------------------------------------------------------

## 5.2 Vaccination delay

Delayed vaccination arms are assigned predefined delays.

The delay represents the operational compromise required to generate a
comparator group.

------------------------------------------------------------------------

## 5.3 Vaccine effectiveness

True biological VE is varied:

  Parameter            Values
  ----------- ---------------
  VE            0.3, 0.5, 0.7

The trial then estimates VE from simulated outcomes.

------------------------------------------------------------------------

# 6. Outcomes

## 6.1 Primary statistical outcome

The principal statistical endpoint is:

\[ VE = 1 - `\frac{risk_{vaccinated}}{risk_{comparison}}`{=tex} \]

The model records:

-   estimated VE;
-   standard deviation of VE estimates;
-   precision metric.

Precision is defined as:

\[ precision = `\frac{1}{1 + SD(VE)}`{=tex} \]

Higher precision indicates more reliable estimation.

------------------------------------------------------------------------

## 6.2 Public health outcomes

The model estimates:

-   direct cases prevented;
-   contacts-of-contacts cases prevented;
-   total cases prevented;
-   vaccine doses used;
-   cases prevented per 1000 doses.

------------------------------------------------------------------------

# 7. Simulation caching architecture

Because many design scenarios share the same epidemic assumptions,
epidemic generation is separated from intervention evaluation.

Cached objects contain:

-   ring structures;
-   latent SAR assignments;
-   epidemic realisations.

They do not contain:

-   allocation assignment;
-   vaccination delay;
-   VE assumptions.

This ensures different trial designs are compared against identical
epidemic conditions.

Cache keys include:

    rings
    simulation number
    SAR
    SAR CV

Example:

    cache_400_2000_0.05_0.25.rds

represents:

-   400 rings;
-   2000 epidemic simulations;
-   SAR 0.05;
-   SAR CV 0.25.

------------------------------------------------------------------------

# 8. Run reproducibility

Each model execution generates a unique run directory:

    runs/
        ring_vax_run_<timestamp>/

            config.rds

            cache/
                epidemic cache files

            checkpoints/
                intermediate design outputs

            outputs/
                operating characteristics
                Pareto frontier

The configuration file records:

-   model parameters;
-   SAR assumptions;
-   SAR CV assumptions;
-   ring numbers;
-   simulation numbers;
-   debug status.

------------------------------------------------------------------------

# 9. Analysis framework

The model supports identification of efficient designs.

A design is evaluated on:

## Benefits

-   cases prevented;
-   efficiency per dose.

## Costs

-   vaccine doses;
-   delay introduced;
-   loss of precision.

Pareto frontiers identify designs where improving one objective requires
sacrificing another.

Importantly, Pareto comparisons are performed within epidemiological
scenarios:

-   same VE;
-   same SAR;
-   same SAR CV;
-   same delay assumptions.

This prevents unrealistic comparisons between fundamentally different
outbreak conditions.

------------------------------------------------------------------------

# 10. Interpretation

The model is intended to answer practical outbreak questions:

1.  How many rings are required to estimate VE with acceptable
    precision?
2.  How much does latent ring-level heterogeneity reduce information?
3.  What is the marginal benefit of moving from balanced randomisation
    toward greater immediate vaccination?
4.  How much additional vaccine is required for extended ring
    strategies?
5.  Under what outbreak conditions is a randomised design operationally
    and ethically justified?

The model therefore evaluates not only whether a design is statistically
efficient, but whether it is a reasonable public health decision.


# 11. Analysis framework and decision outputs

The analysis layer evaluates simulated designs according to two competing objectives:

1. maximising public health benefit;
2. maximising precision of vaccine effectiveness estimation.

Analyses are performed on the completed operating-characteristics output table.

---

## 11.1 Pareto frontier analysis

A Pareto frontier is used to identify designs for which improvement in one objective requires compromise in another.

The primary objectives are:

### Maximise

\[
\text{Total cases prevented}
\]

### Minimise

\[
SD(VE)
\]

A design is considered dominated if another design:

- prevents at least as many cases;
- has equal or lower uncertainty in VE estimation;
- is strictly better in at least one objective.

Pareto comparisons are performed within epidemiological scenarios:

- true biological VE;
- SAR;
- SAR CV;
- delay assumptions.

This prevents inappropriate comparisons between fundamentally different outbreak conditions.

Outputs include:

- Pareto-optimal designs;
- number of Pareto designs per scenario;
- range of cases prevented;
- range of achievable VE precision.

---

## 11.2 Allocation trade-off analysis

The model evaluates the effect of unequal allocation between immediate and delayed vaccination.

For each allocation ratio, the analysis estimates:

- additional cases prevented compared with balanced allocation;
- change in VE uncertainty.

This quantifies the ethical efficiency trade-off:

\[
\text{Additional prevention}
\quad vs \quad
\text{loss of precision}
\]

The objective is not to identify a universally optimal allocation ratio, but to quantify the marginal benefit of increasing immediate vaccination.

---

## 11.3 Vaccination delay analysis

Vaccination delay is evaluated as an operational parameter.

The analysis examines how increasing delay affects:

- total cases prevented;
- VE estimation characteristics.

This reflects the practical decision that delayed vaccination may improve inference but reduces expected prevention.

---

## 11.4 Number of rings required for VE precision

A central analysis estimates the number of rings required to achieve predefined levels of uncertainty.

Precision thresholds are defined using:

\[
SD(VE)
\]

with example thresholds:

- SD(VE) ≤ 0.10;
- SD(VE) ≤ 0.05;
- SD(VE) ≤ 0.02.

For each threshold, the analysis identifies:

- minimum number of rings required;
- achieved SD(VE);
- dependence on:
  - SAR;
  - true VE;
  - allocation ratio;
  - SAR heterogeneity.

This directly addresses the operational question:

> How many rings are needed to estimate vaccine effectiveness with acceptable precision?

---

## 11.5 Latent SAR heterogeneity analysis

The SAR_CV sensitivity analysis evaluates the effect of unobserved variation in transmission risk between rings.

For each SAR_CV scenario the analysis estimates:

- minimum rings required for target precision;
- best achievable precision;
- loss of information due to heterogeneity.

The key hypothesis is:

\[
SAR_{CV} \uparrow
\Rightarrow
SD(VE) \uparrow
\]

because between-ring variation increases outcome variance.

Expected consequences:

- average cases prevented remain approximately unchanged;
- uncertainty in VE increases;
- more rings may be required to achieve the same precision target.

---

## 11.6 Marginal efficiency of additional rings

The analysis calculates the incremental public health benefit of increasing the number of rings:

\[
\frac{\Delta \text{cases prevented}}
{\Delta \text{rings}}
\]

This identifies diminishing returns from expanding the trial.

Outputs describe:

- additional cases prevented per additional ring;
- how marginal benefit changes by:
  - SAR;
  - VE;
  - allocation strategy.

---

## 11.7 Decision tables

Decision tables provide operational summaries under competing priorities.

Two strategies are identified:

### Inference priority

Select designs that achieve near-optimal precision and maximise prevention among those designs.

### Prevention priority

Select designs that achieve near-maximal prevention while minimising uncertainty.

These summaries translate the simulation outputs into practical design choices.

---

## 11.8 Primary analysis outputs

The analysis pipeline produces:

- Pareto frontier tables;
- allocation trade-off figures;
- delay sensitivity figures;
- rings-required tables;
- SAR heterogeneity sensitivity outputs;
- marginal ring efficiency tables;
- operational recommendation tables.

Together these analyses provide a quantitative basis for selecting a ring vaccination study design under outbreak constraints.
