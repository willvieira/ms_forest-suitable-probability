# Introduction

Climate warming poses a significant challenge for several species, particularly for trees whose range shifts struggle to follow temperature warming [@Sittaro2017].
To understand how species will respond to novel environmental conditions arising from climate change, it is imperative to untangle the mechanisms governing their range limits across the environment.
The niche theory predicts that a species will be present in suitable environmental conditions that allow the species to have a positive growth rate [@Hutchinson1957].
From this theory, we can define the geographic distribution of a species as a manifestation of individual demographic rates, such as growth, survival, and recruitment [@holt2009].
By assuming these demographic rates change with the environment, we can predict a species' range limits based on its individuals' performance [@maguire1973niche; @holt2009].

Biotic interaction is undoubtedly an essential driver of demographic rates and, thereby, should affect species range limits.
A recent theoretical framework based on the coexistence theory has been proposed to assess how biotic interactions can scale up to affect range limits [@Godsoe2017].
Formally, this framework evaluates the intrinsic population growth rate when the focal species is rare [@Chesson2000a], both in scenarios where there is no competition (fundamental niche) and when competitive species reach equilibrium (realized niche).
Numerous studies have explored the influence of climate and competition on the distribution of forest trees across their ranges.
For instance, @Ettinger2017 observed in field experiments that neighboring competition constrained individual performance within the range but facilitated better performance outside the range.
Using a dynamic forest model, @Scherrer2020 showed how slow demographic rates and negative competition reduce the uphill migration rate of 16 tree species.
Despite the evidence, the application of this framework to predict the geographic distribution of species based on demographic rates often reveals weak correlations between the performance of tree species and their distribution [@Thuiller2014; @Csergo2017; @bohner2020; @LeSquin2021; @Midolo2021; @Guyennon2023].

One possible explanation for such discrepancy between demographic rates and species distribution is the common practice of assessing performance under average conditions and pointwise estimations, neglecting the associated uncertainty in these estimates.
Variability in demographic rates and subsequently in population growth rate ($\lambda$) arises from two primary sources [@van2020].
The first source is attributed to demographic and environmental stochasticity, wherein individuals exposed to identical conditions may exhibit different responses simply by chance [@Caswell2009].
The second source of variability arises from heterogeneity encountered at various scales.
These differences can manifest between individual stages that motivated the development of structured population models [@Lewis1942; @leslie1945], and extend to large-scale differences in neighboring patches, as described by metapopulation theory [@Levins1969].
A final source of variability in $\lambda$ is the uncertainty associated with the parameters of the demographic models.

Theory predicts that the uncertainty arising from stochastic and heterogenous processes may lead to divergent outcomes in $\lambda$.
Demographic and environmental stochasticity may increase the uncertainty in $\lambda$, consequently increasing the extinction risk, particularly for populations with low performance or low density [@Holt2005; @Gravel2011].
Indeed, demographical stochasticity increased the extinction risk of European forest trees at the hot edge of their distribution [@Guyennon2023].
On the other hand, heterogeneity has been described as a mechanism that can buffer against the stochasticity in demographic rates, thereby increasing population persistence [@milles2023].
This is particularly relevant in nonlinear models, where Jensen's inequality predicts - for convex response functions - that higher levels of demographic and environmental stochasticity increase the average population growth rate [@Koons2009].
Furthermore, demographic and environmental stochasticity influence abundance variation, indirectly impacting $\lambda$ through density-dependence [@May1978; @Terry2022].
A comprehensive understanding of the response of forest trees to climate change requires incorporating the multiple sources of variability arising from spatio-temporal variation and parameter uncertainty.

Here, we use a stochastic Integral Projection Model (IPM) to predict species-level intrinsic population growth ($\lambda$) for 31 eastern North American tree species.
The IPM integrates the growth, survival, and recruitment demographic rates, which vary in response to climate and competition.
By fitting each demographic rate using non-linear hierarchical Bayesian models, we capture parameter uncertainty at both the individual and local population scales.
Additionally, our model naturally accommodates observed spatio-temporal stochasticity in climate and competition.
Then, rather than ignoring these sources of variability, we embrace them into $\lambda$ by defining species performance through a probabilistic framework.
Specifically, we introduce a novel metric called **suitable probability**, derived from the average population growth rate and its associated variability.
This metric determines the probability of a positive population growth rate for a species under specific climate and competition conditions.

Within this new framework, we investigate how a species' suitable probability changes from the center of its distribution toward the cold and hot borders.
We specifically disentangle the relative impacts of climate and competition on changing suitable probability from the center to the borders.
Given that the cold and hot range limits of each species are relative to their distribution across North America, we further examine whether the effect of climate and competition on suitable probability depends on the biogeographic position of the species.
We conclude by discussing a novel theory that uses suitable probability to establish a link between individual demographic rates and metapopulation dynamics.

# Methods

## Population model

We use an Integral Projection Model (IPM) to predict the intrinsic population growth rate ($\lambda$) as a function of climate and competition.
The IPM serves as a mathematical formulation describing the dynamics of a continuous trait distribution ($z$) within a population over time:

$$
n(z', t + 1) = \int_{L}^{U} \, k(z', z, X, \theta)\, n(z, t)\, \mathrm{d}z
$${#eq:ipm}

In our case, the trait $z$ is defined as the tree's diameter at breast height (dbh), constrained between the limits $L$ and $U$.
The continuous distribution $n(\cdot)$ of dbh $z$ of a population at time $t$ transitions to the next time step using a projection kernel ($K$).
The kernel $K$, with covariates $X$ and parameters $\theta$, comprises three demographic submodels:

$$
k(z', z, \theta) = [Growth(z', z, X, \theta) \times Survival(z, X, \theta)] + Recruitment(z, X, \theta)
$${#eq:kernel}

The growth model assesses the probability of an individual of size $z$ at time $t$ transitioning to size $z'$ at time $t+1$.
The survival model determines the probability of an individual with size $z$ at time $t$ surviving to the next time step.
Lastly, the recruitment model determines the number of new individuals entering the population at each time step as a function of total density $z$.

With the defined $K$, we can estimate the intrinsic population growth rate for a determined set of conditions from the covariates $X$ and sampled parameters from the posterior distribution $\theta$.
Specifically, we discretize the continuous kernel $K$ using the mid-point rule [@Ellner2016] and estimate the intrinsic population growth rate using the dominant eigenvalue of the discretized $K$.
This approach is a local approximation of the population growth rate at the initial time steps.

We employed non-linear mechanistic equations grounded in ecological reasoning to formulate the growth, survival, and recruitment demographic models within the IPM.
Each demographic model varies as a function of the mean annual temperature, mean annual precipitation, and stand basal area of larger individuals ($X$).
Each model's parameters ($\theta$) are species-specific, as each model is fitted separately for each species.
Both climate variables influence each demographic model through an unimodal link function, where each model exhibits an optimal climate and niche breadth for temperature and precipitation.
Additionally, density dependence is integrated based on the plot's total basal area of larger individuals.
Stand density affects growth and survival through a linear model, in which two parameters determine the strength of interaction from conspecific and heterospecific competition.
For the recruitment model, the annual ingrowth rate is modulated by conspecific stand basal area, using an unimodal function to account for both the positive effect of seed source and the negative effect of conspecific competition.
Furthermore, the annual survival rate of potential ingrowth individuals decreases linearly with the stand density of heterospecific individuals. 
Finally, the intercept of each growth, survival, and recruitment model incorporates plot-level random effects to control for the variance shared within the plot-year observations.

We use two open inventory datasets from eastern North America: the Forest Inventory and Analysis (FIA) dataset in the United States [@OConnell2007] and the Forest Inventory of Québec [@Naturelles2016].
These inventories, with multiple dbh measurements over time and space, allow us to use the transition information between measurement years for predicting growth, survival, and recruitment rates.
We selected the 31 most abundant species, comprising 9 conifer species and 21 hardwood species, well-dispersed across shade tolerance and succession traits (Supplementary Material 1).
Furthermore, these species are well distributed across the eastern North American gradient.
A detailed description of the data and model development is available in Chapter 2.

## Defining suitable probability

The estimation of $\lambda$ occurs at the local population scale, specifically at the plot level in our study.
Within a given geographic location, such as a specific latitude where several plots are located, the variance of $\lambda$ among those plots arises from spatio-temporal variations in both climate and competition covariates.
For instance, climate stochasticity introduces noise in annual temperature and precipitation, leading to environmental variations. 
Similarly, even with identical climate conditions, two locations can exhibit different community abundance and composition, resulting in variability in the strength of competition. 
Beyond these spatio-temporal environmentally-induced variations, $\lambda$ can still vary due to model uncertainty.

In our study, demographic models track uncertainty at two complementary scales: individual and plot levels.
At the individual level, plots with the same climate and competition conditions may have different $\lambda$ values due to the uncertainty in the demographic models. 
Similarly, even with the same environmental conditions and averaged parameter values (eliminating demographic uncertainty at the individual level), two plots can still yield different $\lambda$ values due to the spatial uncertainty of each demographic model assigned among plots through the random effects.
Therefore, variability in the population growth rate can arise from spatio-temporal variations in both the environment and the parameters.

Given these different sources of variability in $\lambda$, we define the suitable probability as the area under the distribution for $log(\lambda) \geq 0$.
To estimate this, we first determine the cumulative distribution function, $F(x)$, from the generic probability density function, $log(\lambda) = f(t)$, as follows:

$$
F_{\lambda}(x) = P(\lambda \le x) = \int_{-\infty}^{0} f(t)dt
$${#eq:cdf}

This function represents the cumulative distribution from $-\infty$ to $x$.
Subsequently, we define the suitable probability ($\Lambda$) as the inverse of the cumulative distribution function for $x = 0$, assuming $\lambda$ is at the log scale:

$$
  \Lambda = 1 - F_{\lambda}(0)
$${#eq:sp}

## Modeling suitable probability

We can evaluate the suitable probability of a species at various scales, ranging from a single local plot up to several plots in a region.
At the plot level, sources of variability in $\lambda$ stem from parametric uncertainty at the individual level and temporal environmental stochasticity in climate and competition.
When assessing suitable probability at the metapopulation scale by considering multiple plots simultaneously, we can additionally account for spatial variability in the environment and spatial uncertainty in plot-level parameters.

With the exception of parameter uncertainty at the individual level, all other sources of variability exhibit spatial dependence.
This implies that environmental variability (from climate, competition, or both) and parameter uncertainty at the plot level can vary based on their spatial location.
For instance, plots at the borders of the species distribution may experience more temperature variability than those at the center.
Additionally, plot-level parametric uncertainty can be spatially structured, capturing potential features of demographic variability beyond the climatic and competition covariates, such as historical factors or local edaphic conditions. 

At the metapopulation scale, we can model how suitable probability changes across the species' range distribution, considering both fixed climate and competition effects and the underlying spatio-temporal variability.
In our study, we are particularly interested in examining how suitable probability evolves from the center toward the cold and hot ranges. 
For that, we categorized all plot-year observations of a species based on the gradient of mean annual temperature (MAT), divided into cold and hot ranges using the MAT centroid among all plots ($\frac{max(MAT) + min(MAT)}{2}$).
We chose to use MAT instead of latitude because we are interested in the species' climatic niche, although the two variables are highly correlated.
The choice of MAT over latitude is motivated by our interest in understanding the species' climatic niche.
Despite our choice, MAT and latitude are highly correlated for all species.

We assessed suitable probability separately for the cold and hot ranges, employing a linear model to determine the relationship between $\lambda$ and MAT.
The spatio-temporal variability of $\lambda$ arising from environmental stochasticity and parameter uncertainty influences the variance of the linear model.
As this variance may change depending on the range position, we introduce a submodel for the variance of the linear model to be dependent on MAT.
To accommodate potential asymmetry in this variance, we use a Skew Normal Distribution ($\mathcal{SN}$) incorporating an additional parameter ($\alpha$) that can introduce right or left-skewed tails to the variance:

$$
\begin{split}
&log(\lambda) \sim \mathcal{SN}(\xi, \omega, \alpha) \\[2pt]
&\xi = \beta_{1, \xi} \times MAT + \beta_{0, \xi} \\[2pt]
&\omega = e^{\beta_{1, \omega} \times MAT + \beta_{0, \omega}}
\end{split}
$$ {#eq:metamodel}

Here, $\xi$ is the location parameter or the $\lambda$ average, and $\omega$ is the scale representing the variance around the mean.

### Simulations

We computed $\lambda$ for each species based on the plot-year observations in the dataset, considering both environmentally induced variability and parameter uncertainty.
For every observed species-plot-year combination, we incorporated temporal stochasticity in climate conditions by using the mean and standard deviation of mean annual temperature and precipitation calculated from the years between measurements.
For instance, in the case of a plot observed twice, we calculated $\lambda$ for the second observation with climate conditions drawn randomly from a normal distribution with mean and standard deviation defined from climate observations within the year interval.
Similarly, temporal stochasticity in competition arises from variation in abundance and composition between measured years.
By iteratively performing this calculation, drawing parameter values randomly from the posterior distribution, we introduced demographic uncertainty at the individual level.
For each species-plot-year measurement, we replicated the calculation of $\lambda$ 100 times. 
By applying this approach across all plots, we naturally incorporate spatial variation in climate and competition conditions and spatial uncertainty in plot-level parameters.

For each species-plot-year-replication combination, we calculated $\lambda$ under two simulated conditions.
The first condition assumed no competition, with heterospecific competition defined at zero and conspecific competition at a minimal proportion (total population size ($N$) = 0.1).
This condition simulated the population growth rate under the fundamental niche.
The second condition simulated the invasion growth rate, representing the population growth rate when the focal species is rare ($N = 0.1$) and heterospecific competition set to the observed abundance of the competitive species.
This condition simulated the population growth rate under the realized niche.
The code for these computations on High-performance computing is available on GitHub: [https://github.com/willvieira/forest-IPM/tree/master/simulations/lambda_plot](https://github.com/willvieira/forest-IPM/tree/master/simulations/lambda_plot).

After calculating $\lambda$ for each species-plot-year-replication combination, we used these estimates to fit the linear model describing how $\lambda$ and its variability change across the mean annual temperature gradient.
We fitted the species-specific linear model for the hot and cold ranges using the Hamiltonian Monte Carlo (HMC) algorithm via the Stan software [version 2.30.1 @stan2022stan] and the `cmdstandr` R package [version 0.5.3 @cmdstanr].
To reduce the sampling time, we used a sample of 5000 plots for each species to fit the model.
This sample was necessary only for 6 out of the 31 species.
The R and Stan scripts for these computations on High-performance computing are accessible on GitHub: [https://github.com/willvieira/forest-IPM/tree/master/simulations/model_lambdaPlot/](https://github.com/willvieira/forest-IPM/tree/master/simulations/model_lambdaPlot/).

With the fitted models, we leveraged the posterior distribution to estimate the suitable probability of a species for any value of MAT under fundamental or realized niches for the cold and hot ranges.
Specifically, we estimated suitable probability under four different MAT conditions encountered by the species: at the border and the center of each cold and hot range.
We defined the border of the cold range as the minimum observed MAT for the focal species in the dataset, while the hot range was defined as the maximum observed MAT.
The center location is defined as the centroid of MAT for the focal species.
Although the center location has the same MAT for the cold and hot ranges, both are retained because the model is fitted separately for the cold and hot ranges.
Finally, we estimated suitable probability for each location under no competition (fundamental niche) and heterospecific competition (realized niche) conditions, using the empirical cumulative distribution function over 1000 predictive draws.


# Results

### Model fit

We used a linear model to summarize the evolution of the population growth rate ($\lambda$) and its variability across the cold and hot ranges (Equation @eq:metamodel).
Figure @fig:res_example shows the observed $\lambda$ distribution of $\lambda$ and the fit of the underlying model on the mean annual temperature gradient for the species *Abies balsamea*.
Each point represents a plot-year-replication encompassing the complete spatio-temporal sources of variability arising from the stochastic environment and parameter uncertainty.
The black line represents the average fitted model (how $\lambda$ changes with MAT), and the shaded area is the 90th quantile of model uncertainty (how the variability of $\lambda$ changes with MAT).
From this uncertainty, we can deduce the suitable probability.
In this example, in the cold range, the mean and variance of $\lambda$ decrease towards the cold boundary.
By comparing the model under heterospecific competition with that without competition, we can observe a slight decrease in the average at the cold border but a larger decrease in uncertainty (Figure @fig:res_example, bottom left).

![Distribution of stochastic population growth rate ($\lambda$) for *Abies balsamea* over the mean annual temperature gradient for cold (left panels) and hot (right panels) ranges under no competition (fundamental niche) and heterospecific competition (realized niche). The dots represent $\lambda$ over the plot-year-replication combinations. The model's average line and 90% prediction intervals are estimated using 500 draws from the posterior distribution.](https://willvieira.github.io/book_forest-demography-IPM/extinction_risk_files/figure-html/fig-res_example-1.png){#fig:res_example width=100%}

We can estimate the suitable probability using the empirical cumulative distribution approach (Equation @eq:sp) from the linear model predictions.
The Figure @fig:sp-example shows the suitable probability expected over the mean annual temperature of the same species.
As expected, the suitable probability was reduced toward the cold border and was lower under heterospecific competition (yellow curve).
We can also observe that the decrease in suitable probability toward the border is nonlinear, becoming more substantial for heterospecific competition than for the no-competition condition.
The model fit and the estimation of suitable probability across the temperature gradient for all species are presented in Supplementary Material 2.

![Suitable probability of *Abies balsamea* over the mean annual temperature gradient for cold and hot ranges under no competition (green) and heterospecific (yellow). The vertical dotted line represents the range limits of the MAT observed in the dataset.](https://willvieira.github.io/book_forest-demography-IPM/extinction_risk_files/figure-html/fig-sp-example-1.png){#fig:sp-example width=100%}

### Effect of climate and competition on suitable probability

We analyzed the effect of competition on the suitable probability at the border and center of the temperature range distribution for all species.
Figure @fig:sp_comp_vs_nocomp2 compares the suitable probability under no competition to those under heterospecific competition for four locations from the mean annual temperature gradient.
Overall, suitable probability was high among the species, with an average of 0.78.
Among the four locations, species presented a lower suitable probability at the border of the hot range, with an average of 0.67.
Almost all species, in all conditions, are distributed below the identity line (1:1), meaning that heterospecific competition reduces the suitable probability.
In the cold range, temperature decrease toward the border had little effect on the suitable probability for both competitive conditions.
However, at the hot range, however, we can observe a significant shift in suitable probability from the center to the border.
Across the temperature range, there is a monotonic decrease in suitable probability from the cold border toward the hot border.

![Suitable probability for the 31 forest species in their distribution from the cold to the hot border of the mean annual temperature gradient. Note that we omitted the parameter uncertainty of each species in this figure to avoid overlap and increase clarity.](https://willvieira.github.io/book_forest-demography-IPM/extinction_risk_files/figure-html/fig-sp_comp_vs_nocomp2-1.png){#fig:sp_comp_vs_nocomp2 width=100%}

While the negative effect of climate on suitable probability is evident, distinguishing the individual effects of competition remains challenging.
Hence, we further assessed the impact of competition on reducing suitable probability by subtracting the suitable probability under heterospecific competition from the suitable probability with no competition (Figure S1).
Even when isolating the effect of climate, the negative impact of heterospecific competition on suitable probability increases toward the border of the hot range.

### Suitable probability change from center to border

For each species, we estimated the difference in suitable probability between the border and center (for climate) and between heterospecific and non-competition (Figure S2).
We further estimated the relative effect of climate and competition on changing suitable probability from the center to the border of the species distribution (Figure @fig:diff_sp_hist).
Specifically, we quantified the distinct contributions of temperature and competition to changes in suitable probability.
Positive values of the relative difference signify an increase in suitable probability from the center towards the border, negative values indicate a decrease, and zero implies no net change.
In the hot range, most species exhibited negative values of the relative difference in suitable probability, indicating a decline in suitability due to higher temperature and intensified competition toward the hot border of the species distribution.
At the cold range, most species showed positive values for the competition effect, meaning a reduction in the effect of competition toward the cold border.
However, the climate effect at the cold range showed a less clear pattern, with some species experiencing an increase and others a decrease in suitable probability as temperature decreased toward the cold range.
Overall, the relative difference in suitable probability from the center toward the border was more influenced by climate rather than competition.

![Difference in suitable probability for climate and competition effects over the cold and hot ranges. Negative values denote a decrease in species suitable probability from the center towards the distribution border, while positive values indicate an increase. Specifically, a negative value for climate at the hot (or cold) range signifies a reduction in suitable probability as temperature rises (or falls) towards the border.](https://willvieira.github.io/book_forest-demography-IPM/extinction_risk_files/figure-html/fig-diff_sp_hist-1.png){#fig:diff_sp_hist width=100%}

In our investigation of whether suitable probability decreases simultaneously at both borders (indicating a unimodal shape) or exhibits a linear pattern from one border to the other (Figure S3), we found distinct patterns under the climate and competition effects.
For the climate effect, the majority of species demonstrated a decrease in suitable probability at one border while the other remained unchanged.
Additionally, a few species displayed a clear linear pattern of decreasing suitable probability from the cold to the hot border, with only one species (*Betula papyrifera*) showing a pronounced unimodal shape.
Conversely, under the competition effect, most species exhibited a decrease in suitable probability at the hot border and an increase at the cold border, indicating a linear rise in the impact of competition from the cold to the hot border of the distribution.

### Suitable probability across the mean annual temperature gradient

We further investigated the potential dependency of climate and competition effects on the geographical location of the species.
Specifically, we hypothesized that the impact of climate might be more pronounced for species distributed in hotter temperatures than those in colder temperatures.
In Figure @fig:sp_diff_over_MAT, we evaluated the difference in suitable probability for climate and competition effects along the mean annual temperature gradient.
Notably, the only condition showing a significant trend was the climate effect in the cold range (blue group at the top panel).
Not surprisingly, suitable probability decreased as temperature declined toward the cold range, but this pattern was observed only for species located in colder conditions.
Conversely, the suitable probability increased when species were situated in hotter conditions as temperature decreased toward the cold range.

![Relative difference in suitable probability between the center and border for climate and competition for 31 species located over the mean annual temperature gradient. Species points are grouped by a Multivariate Normal Density function with 75% probability.](https://willvieira.github.io/book_forest-demography-IPM/extinction_risk_files/figure-html/fig-sp_diff_over_MAT-1.png){#fig:sp_diff_over_MAT width=100%}

# Discussion

Understanding the mechanisms shaping species distribution is imperative to face ongoing global changes.
We acknowledged and integrated various sources of variability in the population growth rate of forest trees, contributing to an improved understanding of forest dynamics in an uncertain world.
Introducing a novel metric, we quantified the relative impacts of climate and competition on the change in suitable probability across species distributions.
Our findings revealed a nearly linear reduction in suitable probability from the cold to hot borders.
Notably, the predominant influence on the relative difference in suitable probability from the center toward the border was attributed to climate rather than competition.
These results, supported by a novel approach accounting for uncertainty, enhance our understanding of the nuanced interplay between climate and competition across species ranges.

The suitable probability was high across all species and range locations, with only around 5% of all species-location combinations having a suitable probability below the 0.5 threshold.
This is primarily attributed to most species exhibiting a high positive population growth rate across their current range distribution.
Additionally, the spatio-temporal variability in the environment and the parameter uncertainty in the plot may contribute to the elevated average population growth rate due to nonlinear averaging.
This aligns with theoretical [@Schreiber2009] and empirical [@Crone2016] work indicating that spatial heterogeneity leads to a higher population growth rate.

Competition had a significant effect in reducing suitable probability across all range locations, contributing to the ongoing debate surrounding its significance in setting range limits.
Despite several studies emphasizing the effect of competition compared to climate on the demographic rates of forest trees [@Zhang2015;@Kaber2021;@LeSquin2021], debates persist regarding whether this effect at the local scale translates to the biogeographic distribution of species [@Soberon2007;@CopenhaverParry2017].
Our findings support the @Godsoe2017 hypothesis and a growing body of evidence [@Scherrer2020;@Shi2020;@Paquette2021;@Lyu2022] showing that the effect of competition on the intrinsic population growth rate can indeed set range limits.

The decline in suitable probability from the cold to the hot border suggests a predominantly linear, rather than unimodal, performance pattern across the temperature range for most species.
This result is consistent with reduced population growth rates in North American [@schultz2022;@LeSquin2021] and European [@Guyennon2023] forest trees, except for the contrasting pattern observed by @Purves2009.
The higher suitable probability in the cold range compared to the hot range could be attributed to multiple factors.
First, species may still follow their climate niche post the last glaciation, explaining why the current cold range limit does not align with the expected niche distribution [@Svenning2007], potentially leading to a colonization debt [@Talluto2017].
Notably, four of the six species exhibiting a significant decrease in suitable probability from the center toward the cold range were already at the extreme cold observed in the dataset (Figure @fig:sp_diff_over_MAT).

Second, our model might be overlooking crucial drivers of species performance, despite capturing a substantial amount of variation from parameter uncertainty at the plot level.
Factors such as the impact of extreme temperature and precipitation on phenology can influence tree range limits [@Morin2007].
Beyond covariates and plot-level uncertainty, incorporating temporal uncertainty at the plot level, accounting for spatio-temporal covariance, could likely capture additional sources of variation in demographic rates.
While our approach considers temporal stochasticity in climate and competition, which affect species range size [@Holt2022], there remains temporal variation in demographic rates beyond the covariates.
This variability, possibly captured with random effects at the plot level, can influence range limits based on the degree of temporal autocorrelation and its relationship with the range [@Benning2022].
For instance, an empirical study on perennial herbaceous species demonstrated that temporal environmental stochasticity reduced the population growth rate relative to the average [@Crone2016].
In our study, this temporal variability is particularly relevant for survival (due to disturbance) and recruitment (due to phenology) rates because, in addition to having high temporal variability [@Clark1999a; @Leite2023], they represent the most significant drivers of population growth rate (Chapter 2).

The effect of competition, similar to climate, increased from the center towards the border of the hot range, contrary to @Kunstler2021, who found no difference in the competition effect between the center and border of the species.
Additionally, our results deviate from the Species Interactions-Abiotic Stress Hypothesis, predicting a stronger competition effect in less stressful climate conditions [@Louthan2015].
When considering the relative position of the species across the temperature gradient, only the effect of climate at the cold range changed with temperature.
This indicates that most species have a similar or higher suitable probability at the border of the cold range compared to their center distribution.
We further tested whether the species' range size affects the relative difference in suitable probability; while the absolute values change, the pattern among the species remains unchanged.

The climate gradient of temperature had a more significant effect than competition in changing the suitable probability of forest trees.
This means that mean annual temperature, along with all latent variables, better explains how suitable probability changes across the temperature range.
The choice of using only mean annual temperature as an explanatory variable in the metamodel can be improved.
For instance, the model could be built accounting for mean annual temperature and precipitation to predict the complete two-dimensional distribution of the species' climate niche.
Plot random effects could be further used to account for the nestedness of the data design, allowing the proper separation of the total variance of the metamodel into variance arising from individual- and plot-level demographic uncertainty.
While we have assumed climate variability as independent and identically distributed random variables, this assumption can be relaxed to include temporal autocorrelation.
Autocorrelated environmental fluctuation can significantly change a species' range limits due to nonlinear averaging [@Benning2022;@Holt2022].
Lastly, although coexistence theory assumes that the abundance of competitors is at equilibrium [@Chesson2000a], testing this assumption remains practically impossible.

Despite room for improvement in our study, there is a growing body of evidence indicating a mismatch between performance and occurrence [@McGill2012;@Thuiller2014;@Csergo2017;@bohner2020;@LeSquin2021;@Midolo2021;@Guyennon2023].
Our approach can better capture the nuanced effect of climate and competition along with the spatio-temporal variation in $\lambda$, yet it was not enough to fully predict tree range limits.
Since species distribution is influenced by processes at multiple scales  [@McGill2010;@Heffernan2014], it is challenging to rely on a single individual-level performance metric to predict it all [Evans2016].
For instance, dispersion plays a crucial role in changing species distribution at larger spatial scales, either reducing its extent due to limited dispersal or increasing it through source-sink dynamics [@Pulliam2000].
We propose that our novel metric, suitable probability, can be a key unifying factor linking local and landscape scales.

- While forest trees vary in their frequency of occurrence across their distribution gradient, their relative abundance remains constant when present [@Canham2010].
- This suggests that forest distribution should be assessed through colonization and extinction patch dynamics rather than local performance [@Canham2017]. 
- However, instead of attempting to define the best scale for assessing species distribution, we propose using suitable probability to integrate local demographic dynamics with metapopulation theory.
- Describe the metapopulation theory [@Levins1969].
- The extension to include colonization and extinction rates as functions of the environment allows us to assess range limits  [@Holt2005].
- There is an implicit assumption that when a patch is not occupied, it is necessarily available for colonization.
- We relax this assumption and determine the amount of patch availability using the suitable probability.
- When the population growth rate and its variability are high, suitable probability equals 1, and all non-occupied patches are available.
- However, as $\lambda$ decreases or its variability increases, the suitable probability decreases, reducing the proportion of non-occupied patches available for colonization.

# References
