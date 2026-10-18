#import "@preview/lambda:0.1.1": thesis, note-block
#import "@preview/abbr:0.3.1"

// ─────────────────────────────────────────────────────────────────────────────
// Dynamics of Marine Biodiversity and Conservation — Practical report
// Pierre Noire time series (1977-2007), Bay of Morlaix
// Max 15 pages incl. figures; structured as Intro / Methods / Results / Discussion.
// No references required by the assignment, so the bibliography block is left out.
// ─────────────────────────────────────────────────────────────────────────────
#show: thesis.with(
  // ── Document metadata ──
  title: "Twenty Years After: Diversity Dynamics of the Pierre Noire Macrobenthic Community Following the Amoco Cadiz Oil Spill",
  short-title: "Dynamics of Marine Biodiversity and Conservation",
  authors: ("Finn Linus Krauss",),

  // ── Cover page ── (disabled, see toggles below)
  date: "October 2026",

  // ── Front-matter info bar ──
  email: ("finn.linus.krauss@imbrsea.eu",),

  // ── Abstract ──
  abstract: [],
  keywords: "",

  // ── Feature toggles (set any to false to disable) ──
  cover: false,
  show-abstract: false,
  show-outline: true,
  header-footer: true,
  columns: 1,

  // ── Page layout and body text (overrides the template defaults of 1.25 cm side margins and 9 pt) ──
  theme: (margin: (left: 2cm, right: 2cm, top: 2cm, bottom: 2cm), font-size: 11pt),)

// ── Abbreviations ──
#show: abbr.show-rule
#abbr.make(("SR", "species richness"),
  ("Q", "Rao's quadratic entropy"),
  ("FRic", "functional richness"),
  ("FEve", "functional evenness"),
  ("FDiv", "functional divergence"),
  ("FDis", "functional dispersion"),
  ("CWM", "community-weighted mean"),
  ("PCoA", "principal coordinate analysis"),
  ("nMDS", "non-metric multidimensional scaling"),
  ("UPGMA", "unweighted pair group method with arithmetic mean"),
  ("SD", "standard deviation"))

// ── Layout ──
#show heading: it => {
  v(1.6em, weak: true)
it
  v(0.8em, weak: true)}

#show figure: set block(above: 1.4em, below: 1.4em)
#show figure.where(kind: table): it => block(above: 1.4em, below: 1.4em, stack(
  dir: ttb,
  spacing: 0.8em,
  it.caption,
  it.body,))

// ═════════════════════════════════════════════════════════════════════════════
// CONTENT
// ═════════════════════════════════════════════════════════════════════════════

= Introduction

The _Amoco Cadiz_ oil spill of 1978 is an acute example of a direct, large-scale human pressure, and a case study of a community's response and subsequent recovery. It released about 220 000 t of oil, of which between 10 000 and 92 000 t became trapped in subtidal sediments @dauvinFineSandAbra1998. Hydrocarbons were found in the sediments of the Bay of Morlaix two weeks after the wreck.

The Pierre Noire station lies at 17 m depth in the eastern Bay of Morlaix, in a fine sand _Abra alba_ community. At the time, this community was holding the richest benthic amphipod fauna that was recorded in European seas. It also included dense populations of _Ampelisca_ @dauvinFineSandAbra1998. Amphipods lack pelagic larvae, which results in highly isolated populations. In addition, the fine sand communities of the English Channel are separated from each other by other types of sediment, so the population at this site is strongly independent. Sampling at the station began in 1977, just one year before the spill. After the spill, the _Ampelisca_ populations that were highly dominant before almost disappeared and total abundance fell sharply. It took until 1989 for _Ampelisca_ to recolonise the station, and abundance reached pre-spill levels only in 1993.

The original analysis of the series by #cite(<dauvinFineSandAbra1998>, form: "prose") followed abundance, biomass and production from 1977 up to 1996. Today diversity is considered to have several interlinking components, including richness, dominance and functional composition, which can show different facets of the same disturbance and recover at different speeds. The present series now covers 31 years, 11 more than in the original study, and includes a second decline in the last few years. It can therefore be used to follow these components through the spill, the recovery and the period after it.

This report describes the diversity of the Pierre Noire macrobenthic community from 1977 to 2007 with three aims: to describe alpha, functional and beta diversity through time and to identify whether richness or dominance marks the transitions between community phases; to test which diversity indices carry independent information and which overlap with others; and to test whether broader taxonomic resolution (taxonomic sufficiency, @ellisTaxonomicSufficiencyPollution1985) or a single taxonomic group (surrogacy) can replace full species-level identification, to increase the feasibility of long-term monitoring.


= Materials and Methods

== Data and software

The presented analysis used the Pierre Noire macrobenthic time series from the Bay of Morlaix, which comprises 31 annual community samples ranging from 1977 to 2007 @dauvinFineSandAbra1998. The dataset was provided as course material for the course "Dynamics of Biodiversity and Conservation" at the Roscoff Biological Station. It contained a species-by-year abundance matrix (381 species), a taxonomic classification (genus, family, order, class and phylum) and fuzzy-coded biological traits for 224 taxa (8 traits, 36 modalities). Functional indices were calculated from the abundances of the taxa for which trait information was available. Each year was treated as one sampling unit. All analyses were done in R with the packages `vegan`, `FD`, `strucchange` and `ggplot2`, and this manuscript was compiled with Typst.

#note-block(title: "Data availability")[The underlying data, R code, resulting plots and the manuscript's .typ file are available in this project's GitHub repository: #underline(link("https://github.com/finnbambus/Dynamics-of-Marine-Biodiversity-and-Conservation---Report"))]

== Alpha diversity

For each year, the following biodiversity indices were calculated: @SR, the Shannon index ($H'$), the Simpson index, Pielou's evenness and average taxonomic distinctness ($Delta^+$) @clarkeTaxonomicDistinctnessIndex1998. Because $Delta^+$ is based on presence/absence, it is inherently independent of abundance. Hill numbers @hillDiversityEvennessUnifying1973 were used to express richness, Shannon and Simpson on a common scale (Appendix, @fig:hill), but they were not part of the main analysis because the rescaling added little and made the results less intuitive to read. Pairwise Spearman rank correlations, including correlations with total annual abundance, were used to identify redundant indices, meaning indices whose trends were so strongly correlated that they added no information, as well as connections driving each index. Where no p-value is given, it was far below the conventional threshold of 0.05. Of each group of strongly collinear indices, only one representative was investigated further. Dominance structure was described with k-dominance curves @clarkeApproachStatisticalAnalysis and by the total abundance of the five most abundant species of each year.

=== Functional diversity

Functional diversity was calculated from the fuzzy-coded trait table with the `dbFD` function @laliberteDistancebasedFrameworkMeasuring2010, giving equal weight to each trait, in a trait space of five @PCoA axes. @Q, @FRic, @FEve, @FDiv and @FDis were calculated from the same data. All except @FRic:s are abundance-weighted by construction @raoDiversityDissimilarityCoefficients1982 @villegerNewMultidimensionalFunctional2008. Redundancy between indices, and their relationship to total annual abundance, was assessed with Spearman correlations, as for the alpha indices. In addition, @CWM trait values were calculated for body size, life span, feeding mode, development type, movement type and tolerance to disturbance. The five recurrently dominant species (_Ampelisca sarsi_, _A. armoricana_, _A. tenuicornis_, _Spio decoratus_ and _Paradoneis armata_) were excluded when calculating the @CWM:pla, so that these described the remaining community instead of the abundance of the dominant amphipods.

== Beta diversity

Changes in community composition over time were analysed with abundance data, which were $log(x + 1)$-transformed to reduce the influence of the dominant taxa. Bray-Curtis dissimilarities between all pairs of years were visualised using @nMDS in two dimensions, and the points were connected in chronological order to show the community trajectory. The same dissimilarity matrix was clustered hierarchically using @UPGMA, and the tree was cut into five groups, which were interpreted as successive community phases after checking that they formed coherent groups of years in both the tree and the @nMDS. The phases were then used to ask which diversity component marked each transition. Because @SR reflects the number of species and $H'$ combines richness with evenness, a transition matched by a breakpoint in @SR indicates a change in richness, and one matched only by a breakpoint in $H'$ indicates a change in dominance or evenness. Breakpoints in the mean of @SR and of $H'$ were searched with the `breakpoints` function of the `strucchange` package, with a minimum segment length of four years and the number of breakpoints chosen by the Bayesian information criterion. A phase boundary was attributed to a series when the last year of the earlier phase fell inside the 95% confidence interval of one of its breakpoints. Mean richness, $H'$, total abundance and the leading taxa were then compared between phases to describe what changed. A break at the spill could not be tested because the series started only one year before it. To see whether the pattern reflected changes in abundance or in species identity, the ordination was repeated with presence/absence data using the Jaccard and Sørensen coefficients (Appendix, @fig:nmds-jaccard).

== Taxonomic sufficiency and surrogacy

To test taxonomic sufficiency @ellisTaxonomicSufficiencyPollution1985, species abundances were aggregated to genus, family, order, class and phylum. At each level, @SR:s and Shannon were calculated and correlated with their species-level values. To check sufficiency for beta diversity, the Bray-Curtis dissimilarity matrix and @nMDS:s were also recalculated at each taxonomic level. Similarity to the species-level pattern was quantified with a Mantel test and by visual comparison of the ordinations. Taxonomic surrogacy was tested for three groups defined from the given taxonomy: Polychaeta (113 species), Mollusca (93 species) and Amphipoda (86 species). For each group, @SR:s and Shannon were correlated with the whole-community values, and the group's Bray-Curtis matrix was compared with that of the whole community with a Mantel test. A group was considered a useful surrogate if it was significant for all three measures and explained a large share of the variance ($R^2$).

== Statistics

Relationships between indices were tested with Spearman rank correlations because the series were not normally distributed. Linear regressions were used for the trend of $Delta^+$ over time and for the sufficiency and surrogacy comparisons, and all Mantel tests used 999 permutations.

= Results

== Alpha diversity

Shannon, Simpson and Pielou were near-collinear ($rho^2 = 0.87$-$0.94$, all $p < 0.001$; @fig:alpha f), so Simpson and Pielou added no information beyond Shannon and were not interpreted further. @SR was independent of all three (Shannon $rho^2 = 0.02$, Simpson $0.07$, Pielou $0.11$; all n.s.), and $Delta^+$ was uncorrelated with all four other indices (n.s.). The results therefore focus on @SR, Shannon and $Delta^+$, which captured three independent axes of variation.

Species richness ranged from 56 to 109 species per year (with a mean of 78 and a @SD of 13) and showed no clear visual long-term trend (@fig:alpha a). There was also only a weak spill signature right after the disruption: 68 species were recorded in 1977, 59 and 57 in 1978 and 1979 (a variation that falls within @SD), and then a jump to 81 in 1980. The highest @SR (109) occurred in 1991, and lowest values were found again in 2002 (56) and in 2005-2007 (58-67). The time series stopped during a downward trend, so the observed local minimum may have been surpassed in the following years. While the correlation between @SR and total abundance was highly significant, it had only weak explanatory power ($rho^2 = 0.20$, $p = 0.012$; @fig:alpha f).

#figure(
  image("figures/alpha_diversity.png", width: 100%),
  caption: [Alpha diversity of the Pierre Noire macrobenthos by sampling year: panels a and b are the indices used in the main analysis:
    (a) species richness (SR), 
    (b) Shannon index $H'$, 
    (c) Simpson index $1 - D$,
    (d) Pielou evenness $J$, 
    (e) average taxonomic distinctness $Delta^+$. 
    (f) Spearman correlation matrix between SR, Shannon, Simpson, Pielou, $Delta^+$ and total abundance (Abund.); numbers are $rho^2$, colour shows significant pairs (darker meaning a correlation closer to ±1), and grey fields indicate non-significance ($p gt.eq 0.05$).],
) <fig:alpha>

#figure(
  image("figures/dominance_structure.png", width: 100%),
  caption: [Dominance structure of the Pierre Noire macrobenthos.
  (a) k-dominance curves (cumulative relative abundance by species rank, log rank on the x-axis) coloured by year. 
  (b) Stacked abundance of the five most abundant species of each year (coloured bars: the five recurrently dominant species and "Other", i.e. species that were in the top five in a given year but not generally). The black line and points show the total annual abundance, overlaid on the most abundant species. The gap between bars and line shows the combined
    contribution of all remaining species. Abundance in units of
    10 000 individuals.],
) <fig:dom>

Shannon diversity fell into three visible regimes (@fig:alpha b). It was high throughout the 1980s, surging right after the spill in 1978, then dropped in the 1990s to a low plateau between 1993 and 2004, and then rose sharply to the highest values of the series in 2005-2007. The pre-spill year 1977 had the lowest value of the whole series (1.48) and $H'$ jumped to 2.64 in 1978, so the spill was accompanied by a direct increase in $H'$. The visible negative relation between Shannon diversity and total abundance was also confirmed by a highly significant negative correlation ($rho^2 = 0.69$; @fig:alpha f): low abundance coincided with high $H'$ (as visible in 1978-1981 and 2005-2007 in @fig:dom). It also connects directly to dominance, as lows in $H'$ coincide with times of high dominance. In 1993, for example, _Ampelisca sarsi_ alone made up 54% of all individuals, and in 1977 two _Ampelisca_ species made up 84%. This relation is quantified by the strongly positive correlation of Shannon with Pielou's evenness ($rho^2 = 0.94$; @fig:alpha f).

=== Taxonomic diversity

Average taxonomic distinctness varied only little (between 82.0-86.8, mean 85.4) and showed no trend over time (linear regression $R^2 = 0.001$; @fig:alpha e). It also showed no obvious spill signature (85.5 in 1977, 84.6 in 1978). The lowest values occurred in 2002 (82.0) and 1984 (83.5); 2002 was also the year of minimum @SR (56). In the correlation tests, $Delta^+$ was also independent of all other indices: @SR, Shannon (treated as representative for Simpson and Pielou), and total abundance (@fig:alpha f).

=== Abundance and dominance

The k-dominance curves ranged from years in which the five most abundant species already accounted for almost all individuals (84.6-95.4% in 1993-2003) to years that needed more than 20 species to reach 95% of all individuals (most years from 1978 to 1991, and 2005-2007; @fig:dom a). On average, the five most abundant species of each year made up three quarters of all individuals (76%, range 47.5-95.4%), and their summed abundance tracked total annual abundance almost perfectly ($rho^2 = 0.98$; @fig:dom b). These were three _Ampelisca_ species (_A. sarsi_, _A. armoricana_, _A. tenuicornis_), as well as _Spio decoratus_ and _Paradoneis armata_.

=== Functional diversity
#figure(
  image("figures/functional_diversity.png", width: 100%),
  caption: [Functional diversity of the Pierre Noire macrobenthos by sampling year: panels a-c are the indices used in the main analysis: 
  (a) Rao's quadratic entropy, 
  (b) functional richness (FRic), 
  (c) functional evenness (FEve), 
  (d) functional divergence (FDiv), 
  (e) functional dispersion (FDis). 
  (f) Spearman correlation matrix
    among the five indices and total abundance (Abund.); numbers are $rho^2$, colour shows significant pairs (darker meaning a correlation closer to ±1), and grey fields indicate non-significance ($p gt.eq 0.05$).],
) <fig:fd>

@FDis was practically identical to @Q ($rho^2 = 0.99$), and @FDiv also strongly mirrored both ($rho^2 ≈ 0.85$; @fig:fd f). The results therefore focus on @Q, @FRic and @FEve.

Rao's Q ranged from 0.70 in 1994 to 3.61 in 2007 and strongly followed the pattern of Shannon ($rho^2 = 0.91$; @fig:fd a). It was high in the 1980s, fell to a low plateau between 1993 and 2004 and rose sharply after 2004. As in Shannon, the pre-spill year was low and the following year high. As expected from its close connection to Shannon, Q was unrelated to @SR ($rho^2 = 0.03$, n.s.) and negatively related to total abundance ($rho^2 = 0.71$).

Functional richness ranged from 4.22 to 8.50 (@fig:fd b). Its lowest value occurred in 1978, when it dropped by ~30% from the year before; a second low was reached in 2002 (4.99) and the highest values in 1999-2000. Unlike $Delta^+$, @FRic therefore showed a clear minimum at the time of the spill. It was independent of @Q, Shannon and $Delta^+$ (all n.s.), and instead correlated with @SR ($rho^2 = 0.67$) and, weakly, with total abundance ($rho^2 = 0.20$, $p = 0.012$).

Functional evenness varied within a narrow range (between 0.42 and 0.55), with the highest values in 2005 and 2007 (@fig:fd c). It was not correlated with @Q, Shannon, $Delta^+$ or Pielou, but showed a weak negative correlation with @FRic ($rho^2 = 0.24$, $p = 0.005$) and with total abundance ($rho^2 = 0.26$, $p = 0.003$).

#figure(
  image("figures/cwm_traits.png", width: 100%),
  caption: [Community-weighted means (CWM) of six selected functional trait classes per sampling year (size, life span, feeding mode, development type, movement type, tolerance to disturbance), calculated from fuzzy-coded traits, excluding the five most dominant species. Each line is the share of the community carrying that trait category.],
) <fig:cwm>

Community-weighted means (@fig:cwm; dominant species excluded) were variable for most traits. In the mid-1980s, a single trait class peaked in five of the six traits at the same time: very small body size, short life span, surface deposit-feeding, direct development and sensitive taxa all reached their maximum in 1984, a period of low total abundance. Burrowers, the dominant movement type, peaked earlier, in 1980. For tolerance to disturbance, the share of sensitive taxa was high in the mid-1980s, fell in the early 1990s to about half of its earlier level, and rose again from 2000 to a peak in 2005, before decreasing again in 2006-2007. Feeding mode shifted from a dominance of surface deposit-feeders in the 1980s towards a more even share of surface deposit-feeders and predators/scavengers from the 1990s onwards. This reverses the crash in the share of predators/scavengers during the spill, which was followed by the takeover by surface deposit-feeders described above. Development type and life span showed similar patterns to each other. Short-lived taxa dominated in the 1980s, whereas medium-lived taxa were more frequent in the 2000s. Direct development was most frequent in the mid-1980s. Burrowers were the dominant movement type throughout, while the share of tube dwellers peaked around 1990 and 2000 and was low in 2005-2007. Body size shows an increase in very small organisms in the years following the spill, a trend that normalises again after a couple of years of dominance.

== Beta diversity

The @nMDS of Bray-Curtis dissimilarities (stress = 0.156) showed the community moving along a directional trajectory through time (@fig:beta a). Cutting the @UPGMA tree into five groups gave coherent, naturally chronological periods that were also separated in the ordination (@fig:beta a, c): a pre-spill phase (1977), a post-spill crash (1978-1981), a post-spill remnant phase (1982-1992), recovery (1993-2004) and a millennium decline (2005-2007). Total annual abundance fell from 8187 individuals in 1977 to between 1758 and 4505 during the crash, was highest in the recovery phase (up to 21173) and dropped to 1022-1664 in 2005-2007. On presence/absence data, the Bray-Curtis trajectory was lost (Jaccard and Sørensen gave identical configurations, stress = 0.224; @fig:nmds-jaccard), indicating that the phases reflected the abundance of taxa rather than a change in species identity.

To investigate this further, a breakpoint analysis was used to find out whether each transition between these phases was driven by changes in composition or in @SR. It found two breaks in @SR (after 1980 and after 2001) and two in $H'$ (after 1989 and after 2003; @fig:beta b). Each phase change in the Bray-Curtis dissimilarity matrix matched at least one of them. The break between the post-spill crash and the post-spill remnant is present as an @SR break (confidence interval 1977-1984) but not in $H'$. Mean richness rose from 69 to 86 species, while $H'$ stayed stable. Next, the break between the post-spill remnant and recovery matched only $H'$, with a break between 1988 and 1992. Here $H'$ dropped from 2.47 to 1.73 while richness stayed stable. During that time _Ampelisca sarsi_, _A. armoricana_ and _A. tenuicornis_ together made up 83% of individuals. Finally, the last break matched breakpoints in both indices (@SR between 1998-2005 and $H'$ between 1998-2004). While the change was present in both indices, the trends of the following phase were opposing: @SR fell to 61 species while $H'$ rose to 2.96.

#figure(
  image("figures/beta_diversity.png", width: 100%),
  caption: [Beta diversity of the Pierre Noire macrobenthos.
    (a) nMDS ordination of annual community composition (Bray-Curtis dissimilarity, log(x + 1)-transformed abundances), years are connected in chronological order, labels mark the first year of each phase, colours show the five UPGMA-derived phases. 
    (b) Phases with their years, range of total annual abundance (individuals) and the series (@SR, Shannon or both) with a breakpoint whose 95% interval includes the years of the phase change (n/a: spill not testable). 
    (c) UPGMA dendrogram (Bray-Curtis), branches and labels coloured by the same five phases as a).],
) <fig:beta>

== Taxonomic sufficiency and surrogacy

=== Taxonomic sufficiency

The family level was the coarsest level that retained most of the information present at the species level (@fig:suff). Shannon was covered well at all levels, even phylum ($r$ $gt.eq$ 0.86). @SR was well represented only down to the family level ($R^2$ = 0.88), while the order and class levels reached only $R^2$ ≈ 0.46 and the phylum level $R^2$ = 0.33. The ordination pattern followed the same gradient: the Mantel $R^2$ between species-level and higher-level Bray-Curtis matrices decreased from 0.90 at genus and 0.84 at family level to 0.58 at order, 0.36 at class and 0.41 at phylum level, a trend visible in the representation of the original chronological pattern. With 176 families instead of 381 species (46% of the original taxa number), the family level still explained most of the variance in @SR (88%), Shannon (96%) and community dissimilarity (84%).

#figure(
  image("figures/sufficiency.png", width: 100%),
  caption: [Taxonomic sufficiency. 
  (a-d) Species richness (SR) at different taxonomic levels against species-level SR (line represents a linear model with a 95% confidence band). 
  (e) Regression of SR at each taxonomic level to the species-level SR ($r$ slope of the correlation, $R^2$ explained variance). 
  (f) Same as e) for the Shannon index. 
  (g-k) Bray-Curtis nMDS ordination (log(x + 1)-transformed data) at different taxonomic levels, years connected in chronological order and coloured by year phase (see @fig:beta). 
  (l) nMDS stress at each taxonomic level (stress at species level, 0.156) and the variance in species-level Bray-Curtis dissimilarity explained by each level ($R^2$, squared Mantel correlation with 999 permutations, all highly significant with $p = 0.001$, the smallest value possible with 999 permutations).],
) <fig:suff>

=== Taxonomic surrogacy

Of the three tested candidate groups, only Amphipoda emerged as a significant surrogate, with significant p-values for all three indices (@tab:surrogacy). It could explain 62% of the variance in @SR and 68% in community dissimilarity, but only 40% in Shannon. Polychaeta did better on two of the indices (75% for @SR and 73% for Bray-Curtis), but did not track Shannon (n.s.) at all. Mollusca was weaker on all measures and also not significant for Shannon. The @nMDS of the amphipods alone (stress = 0.104) reproduced the overall trajectory of the ordination of all species (@fig:surrogacy c, d).

#figure(
  {
    set text(size: 8.5pt)
    let ns(x) = text(fill: luma(160), x)                       // non-significant rows (p ≥ 0.05)
    table(
      columns: (1.3fr, 0.8fr, 1.4fr, 0.8fr, 0.8fr),
      align: (left + horizon, right + horizon, left, right, right),
      inset: (x: 5pt, y: 2.5pt),
      stroke: none,
      table.header[Group][Species][Index][r][R²],
      table.hline(stroke: 0.5pt),
      table.cell(rowspan: 3)[Polychaeta], table.cell(rowspan: 3)[113], [SR], [0.867], [0.752],
      ns[Shannon], ns[0.080], ns[0.006],
      [Bray-Curtis], [0.857], [0.734],
      table.hline(stroke: 0.5pt),
      table.cell(rowspan: 3)[Mollusca], table.cell(rowspan: 3)[93], [SR], [0.646], [0.417],
      ns[Shannon], ns[0.175], ns[0.031],
      [Bray-Curtis], [0.420], [0.177],
      table.hline(stroke: 0.5pt),
      table.cell(rowspan: 3)[*Amphipoda*], table.cell(rowspan: 3)[86], [SR], [0.785], [0.617],
      [Shannon], [0.629], [0.395],
      [Bray-Curtis], [0.822], [0.676],
      table.hline(stroke: 0.8pt),
    )
  },
  caption: [Surrogacy of three taxonomic groups for the whole community. $r$ is the correlation slope between the group and whole-community values, $R^2$ the share of variance explained. For Bray-Curtis, $r$ represents the Mantel correlation between the two dissimilarity matrices. Rows in grey are not significant ($p gt.eq 0.05$).],
) <tab:surrogacy>

#figure(
  image("figures/surrogacy.png", width: 100%),
  caption: [The best surrogate, Amphipoda, against the whole community for
    (a) species richness and 
    (b) Shannon index (points coloured by year, lines are linear models with 95% confidence bands). 
    nMDS ordinations (Bray-Curtis, log(x + 1) transformed data, with years connected in chronological order) of
    (c) all species and 
    (d) species of the amphipod group only.],
) <fig:surrogacy>

= Discussion

#heading(level: 3, numbering: none, outlined: false)[The spill response]

The _Amoco Cadiz_ spill is the only clearly documented event in the present time series. #cite(<dauvinFineSandAbra1998>, form: "prose") reported that the dominant _Ampelisca_ populations, which are very sensitive to oil, largely disappeared after the spill, together with a large part of the community's biomass and production. The present data supports this in several respects. In the years after the spill, total abundance is low and the community is far less dominated than in the single year before the spill (1977), when two _Ampelisca_ species made up most of the individuals. Shannon diversity ($H'$) and Rao's Q rise sharply during the same years. This rise is mainly the result of the lost dominance: the community becomes more even without becoming richer. A high $H'$ after a disturbance can therefore not be read as a real sign of better ecological condition, or higher diversity.

In the following years, species richness dips only slightly, staying within its normal year-to-year range. Most of the species pool remains at the site after the disturbance, and the spill acts mainly on how individuals are distributed among species. The polychaetes _Spio decoratus_, _Paradoneis armata_ and _Chaetozone_ spp. lead the post-spill crash composition, which matches Dauvin's description of Spionidae and Paraonidae taking over the space previously occupied by amphipods. Functional richness reaches its lowest value in 1978 while @SR barely changes. As @FRic measures the volume of trait space occupied by the species present, its drop suggests that the species missing in 1978 included some at the edges of the trait space.

Unfortunately, all of these deductions rest on only a single pre-spill year. The 1977 sample cannot show whether its strong dominance and low $H'$ were typical for the site or an extreme low or high of a different normal state. The dominated, high-abundance community of the 1990s gives a rough picture of what a pre-spill community could have looked like. While Dauvin found that _Ampelisca_ densities returned to values similar to those before the spill, its comparison also rests on the same single year and supports only basic deductions.

#heading(level: 3, numbering: none, outlined: false)[Richness and dominance as the two drivers of diversity]

The spill response already shows two separate components of diversity: how many species are present (richness) and how evenly the individuals are shared among them (evenness). The biodiversity indices used all fall somewhere along these two axes (@fig:alpha f, @fig:fd f). $H'$, Simpson and Pielou all show strong linear correlation, representing the evenness side of that gradient, while species richness stands alone, representing the richness side. The strong correlation that merges all evenness indices at one end of the gradient, instead of spreading them along it, results from the data being driven by one of the two ends at a time. Among the functional indices, Rao's Q, @FDis and @FDiv form a second tight group that follows the evenness indicators. All of these indices are weighted by abundance, and the near-identity of Rao's Q and @FDis is expected from the way both are constructed.

This group of indices is also strongly tied to total abundance, as the years with high abundance are also the years in which a few _Ampelisca_ species dominate. Here the five most abundant species hold three quarters of all individuals on average, and their summed abundance tracks total abundance almost perfectly (@fig:dom b). Total abundance, dominance and low evenness are therefore close to a single signal. A year with low $H'$ mainly marks a year in which the dominant amphipods are abundant.

Species richness and @FRic form the opposite side of the gradient. They are correlated with each other, which is expected as the volume of trait space grows with the number of species. Both are only weakly related to total abundance and unrelated to the dominance indices. @FEve is the only indicator that is at most weakly related to either of the two sides. Average taxonomic distinctness showed no reliable trend over time and no response to the spill, so it is not discussed further. Because the two axes are independent, @SR or $H'$ alone shows only half of the picture. Richness can change without a change in dominance and the other way around.

#heading(level: 3, numbering: none, outlined: false)[What drives each phase]

Five phases were identified, based on the ordination and clustering of log-transformed abundances. On presence/absence data the chronological trajectory mostly disappears, although a trace of it remains (@fig:nmds-jaccard). This leads to the assumption that the phases are mainly driven by changes in abundance and evenness while the number of species stays largely stable. The breakpoints found in the time series of @SR and $H'$ can give additional insight into which of the two drivers causes each individual transition from one phase to the next, instead of assuming that one of them explains the whole extent of the data (@fig:beta b).

The transition from the post-spill crash to the post-spill remnant phase (1981/1982) is mainly a change in richness. @SR has a breakpoint here while $H'$ trends stay stable. This implies that while some species return to the site, it happens without a matching shift in dominance. The timing roughly fits the recovery of the sediment reported by #cite(<dauvinFineSandAbra1998>, form: "prose"): hydrocarbon concentrations fell to low levels after 1981 and benthic recruitment in 1981 returned to normal conditions. As a result the post-spill remnant phase is the most species-rich period of the series, while (and maybe because) it is not dominated yet. While _Ampelisca sarsi_ is among the most common species, total abundance stays at moderate levels and individuals stay spread across multiple species.

The transition into the recovery phase (1992/1993) is a dominance change. Only $H'$ has a breakpoint here, and richness changes little while three _Ampelisca_ species take over most of the abundance. The return of the amphipods does not appear to displace other species; instead it drives a large surge in total abundance of the community. Dauvin explained the delayed return of _Ampelisca_ with the biology of amphipods, which have no pelagic larvae, low fecundity and poor dispersal, so a decimated population can only recover through recolonisation from other unaffected populations. He dated the recolonisation to after 1989, which agrees with the present $H'$ breakpoint after 1989. The cluster boundary of dissimilarity follows about three years later.

The transition into the millennium decline (2004/2005) involves both drivers, and they move in opposite directions: richness falls while $H'$ rises. This transition is driven by both a loss of species and a loss of abundance of the most dominant species (resulting in increasing evenness).

Each transition has its own signature. A single metric, such as total abundance or $H'$, would fail to show the return of species in the early 1980s, the return of dominance (assuming pre-spill levels) in the early 1990s, and the combined change to decline after 2004 as separate events, each driven by a different combination of drivers.

#heading(level: 3, numbering: none, outlined: false)[Recovery as a return to high dominance]

In the recovery phase the community returns to a dominated, high-abundance state with the lowest $H'$ values since 1977. The amphipod populations recover, while the evenness-weighted indices fall. #cite(<dauvinFineSandAbra1998>, form: "prose") interpreted this dominated state as a stable peaking community, favoured by stable beneficial conditions in the Bay of Morlaix. It reported that _Ampelisca_ reached densities similar to those before the spill within about 15 years. But whether this is a full recovery to a previously present state cannot be decided from the present data. A single pre-spill year can't show what the undisturbed state looked like, nor how much it normally varied. Judging the actual recovery would need a reference state built from several undisturbed years.

#heading(level: 3, numbering: none, outlined: false)[The millennium decline]

The last recorded phase mirrors that of the recovery. During the recovery the dominance of the _Ampelisca_ species builds up and total abundance follows; during the millennium decline both fall to levels even lower than at the start of the recovery. Total abundance falls to the lowest values of the series and evenness rises, with the highest values of $H'$, Rao's Q and @FEve in the series. Richness also falls. The community moves back along the dominance axis that the recovery built up, but unlike after the spill, this time it also loses @SR. In the ordination and the dendrogram the years 2005-2007 do not return to the crash phase, instead they form a new group with highly distinct years.

The decline shares the low abundance and high evenness of the post-spill crash, but the trait composition points to a different process. After the spill, the share of sensitive taxa in the non-dominant community falls and tolerant taxa increase. In 2005 the share of sensitive taxa reaches its highest value since the mid-1980s. The shift towards tolerant taxa that followed the oil spill is therefore missing, which argues against a similar pollution event but does not exclude other disturbances.

The direct causes of the decline cannot be identified with this dataset alone and remain a question for further research. Possible explanations include changes in temperature and hydrography, fishing pressure, and the arrival of non-indigenous species. Testing those requires environmental data in addition to the biological series as well as a longer record. The series ends while the decline is still under way, so it remains open whether the community settles into a new state or returns to amphipod dominance as it did before.

#heading(level: 3, numbering: none, outlined: false)[Functional composition]

The @CWM:pla describe the trait composition of the community without its five dominant species. They add a view that is independent of the amphipod signal. The clearest pattern is in the post-spill remnant phase. In the mid-1980s a single class peaks in five of the six traits at the same time: very small, short-lived, surface deposit-feeding, directly developing and sensitive taxa.

This profile is mixed and covers a range of organisms. Small, short-lived surface deposit-feeders fit the opportunistic polychaetes that #cite(<dauvinFineSandAbra1998>, form: "prose") described as temporarily dominant after the spill. Together with the richness breakpoint, the peak in sensitive taxa suggests that the community of the mid-1980s combines fast colonisers with sensitive species that return to the site after the spill. The peak may also include a release from dominance.

The feeding mode shows the longest-lasting change. Predators and scavengers make up the largest share in 1977 and decline after the spill. Surface deposit-feeders replace them through the 1980s, and from the 1990s onwards both groups contribute about equally. Life span also shifts in a similar way, from short-lived taxa in the 1980s to medium-lived taxa in the 2000s. This fits a community that moves from fast colonisers towards more persistent species as time since the spill increases.

#heading(level: 3, numbering: none, outlined: false)[Taxonomic sufficiency]

Family-level identification appears to retain most of the information on richness, $H'$ and community structure using less than half of the original number of taxa. It could therefore be a cost-effective alternative to species-level identification for this community. Unlike @SR, $H'$ is recovered even at the phylum level, most likely because a few very abundant taxa dominate it, and the dominance of the _Ampelisca_ amphipods stays visible even at highly degraded resolution. Richness is more spread across many rare species, which are merged once families are combined into orders, resulting in @SR and community structure losing much of their information above family level.

Of the groups tested as surrogates, only Amphipoda track richness, $H'$ and community structure together significantly, however the explained variation of $H'$ is far below useful levels, so none of the tested groups is a suitable surrogate. Amphipoda include the dominant _Ampelisca_ species and therefore partly reproduce the dominance signal; in phases of low abundance of these species the representative power could be even lower. Polychaeta and Mollusca only track richness and community structure but not $H'$.

#heading(level: 3, numbering: none, outlined: false)[Limitations]

The data were derived from samples taken at only one site, sampled once per year without replicates, which makes it impossible to separate spatial and within-year variation from the observed temporal trends.

The five phases come from a cut of the @UPGMA tree that was determined visually. The breakpoint tests assumed abrupt shifts with a minimum segment length of four years between breakpoints, so short events and gradual change are not captured. Further, the functional indices rest on trait information for 224 taxa, excluding 157 of the 381 species in the abundance data. No environmental data that could indicate the cause of the millennium decline were analysed, and longer records from before the spill and after the last decline would add further information.

#heading(level: 3, numbering: none, outlined: false)[Conclusion]

Over the 31 years of the series, the Pierre Noire community moved through five phases. The transitions between them are marked by changes in richness, in dominance or in both. The spill caused a loss of dominance first and of species second, the return of the _Ampelisca_ amphipods restored dominance while richness stayed stable, and the late decline removed dominance and species together. Because the two drivers change independently, monitoring needs to observe both; here family-level identification captures most of the signal obtained at species level.

#pagebreak()

#bibliography("Dynamics of Biodiversity.bib", style: "apa")

// ── Appendix ──────────────────────────────────────────────────────────────────
#colbreak()
#set figure(numbering: n => "S" + str(n))
#counter(figure).update(0)
#counter(table).update(0)
#set heading(numbering: none)

= Appendix

== Supplementary Information

#figure(
  image("figures/nmds_presence_absence.png", width: 100%),
  caption: [nMDS ordination of annual community composition on
    presence/absence data: (a) Jaccard and (b) Sørensen coefficients (stress = 0.224 for both), years connected in chronological order. The two
    coefficients give the same configuration. See Figure
    #ref(<fig:beta>, supplement: none)a for comparison with the
    abundance-based (Bray-Curtis) ordination.],
) <fig:nmds-jaccard>

#figure(
  image("figures/hill_numbers.png", width: 80%),
  caption: [Hill numbers of order $q = 0$ (species richness), $q = 1$
    ($exp(H')$) and $q = 2$ (inverse Simpson, $1 slash sum p_i^2$) by
    sampling year, on a common axis (effective number of species).],
) <fig:hill>

== Aknowledgement for the use of gAI

I aknowledge that generative AI tools were used in this project to adapt the R code provided as course material to the present analysis, in the drafting and to outline the manuscript, to proofread it, and to improve translation into English. The analysis design, the interpretation of the results and the conclusions are my work and responsibility, and before including it I checked and edited all AI-assisted content and code.
