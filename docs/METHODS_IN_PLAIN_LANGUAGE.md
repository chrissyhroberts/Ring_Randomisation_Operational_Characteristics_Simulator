# Methods in plain language

## The question

How many rings should a trial recruit, how should those rings be allocated, and how long should delayed vaccination be postponed when the future outbreak is uncertain?

## Why a single sample-size calculation is insufficient

A conventional calculation usually assumes one attack rate, one vaccine effect and one amount of clustering. Those values are not known here. A design that works under a favourable attack rate may produce too few events when transmission is lower or infections are concentrated within a small number of rings.

## The simulation approach

The model creates many hypothetical outbreaks. Each hypothetical outbreak has a specified:

- true vaccine effectiveness;
- average secondary attack rate; and
- intracluster correlation coefficient.

The ICC controls how similar outcomes are within a ring. Ring-specific baseline risks are generated from a beta distribution. Conditional on that ring-specific risk, participant outcomes are Bernoulli. This produces a beta-binomial form of within-ring clustering.

The same operational designs are tested in every epidemiological setting. Each design is repeated many times, allowing the model to measure:

- the average VE estimate;
- the standard deviation of that estimate;
- bias;
- the fraction of replicates returning a valid estimate;
- cases prevented; and
- vaccine doses used.

## Grid adequacy

Each VE–SAR–ICC combination is treated as a possible world. A design is called inferentially adequate in a world when:

- SD(VE) is at most 0.15;
- absolute VE bias is at most 0.05; and
- at least 95% of replicates return a valid estimate.

Grid adequacy is the percentage of possible worlds meeting all three conditions.

All grid points are weighted equally. The resulting percentage is therefore a transparent sensitivity score, not the probability that the future trial will succeed.

## Precision versus prevention

Delayed vaccination preserves a comparison period between the immediate and delayed arms, which helps estimate VE. Earlier vaccination can prevent more disease but shortens that comparison. The model therefore reports a frontier rather than declaring that one design is optimal without reference to resources or policy priorities.

## CoC sensitivity

The CoC risk multiplier specifies how much infection risk a contact-of-contact has relative to a direct contact. This assumption changes the potential public-health benefit and vaccine demand. It does not remove the low-event problem for the primary direct-contact VE estimator.

