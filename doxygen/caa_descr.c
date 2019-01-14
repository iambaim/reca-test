/**
   @namespace caa - Catch At Age
  
*/

/**
  @file
  Doxygen main page.
  Contains the text on the main page including the coding standard.
*/

/**   
   @mainpage ECA - Estimating catch at age

   \author Hanne Rognebakke \htmlonly <a href="mailto:hanne.rognebakke@nr.no">Hanne Rognebakke@nr.no</a>\endhtmlonly
   \author Geir Storvik 

   @section intro Introduction

   %caa is a set of c-routines for MCMC simulations of catch at age. The main routines are to 
   be called from R or Splus but have main programs for testing as well.
   The software is written in C and has been developed by the Norwegian Computing Center (NR)
   as part of seveal project for HI.

   The catch at age is defined as the number of fish caught in each age
   category. For most hauls, only total catch weight is recorded. 
   The standard procedure for estimating the catch at age from this is to first
   find an estimate of the number of fish by dividing the total weight
   by the mean weight. The number of fish in each age category is then
   estimated by multiplying the estimate of the number of fish with an 
   estimate of the proportion of fish in each age category.

   Both the mean weight and the proportion at age will change with
   year, season, gear and regions, but can also change from one haul to
   another. Information of these quantities can be obtained from samples
   where both age and weight, but also length, are observed.
   In \ref Hirst_etal_2004 and \ref Hirst_etal_2005 a hierarchical
   model describing the relationship between age, length and weight is proposed.
   Including length in the model is necessary because some samples only 
   contain length measurements.
   Putting the model into a Bayesian framework, uncertainties in the parameters
   are taken into account in the catch at age estimates.
   CAA is based on this model and framework.
   
   In section \ref sec_estcaa we describe the procedure for estimating catch at age in more 
   details.
   This model for the relationship between age, length and weight
   is described in section \ref sec_mod. 
   
   @section sec_estcaa Estimation of catch-at-age

   Let \f$W_h\f$ be the total weight at haul \f$h\f$. We want to estimate the
   total number of fish in each age groups. Now, \f$T_h\f$, the number of fish in
   haul \f$h\f$, is unknown, but is given by
   \f[
   T_h=\frac{W_h}{\overline{w}_{h,\cdot}}
   \f]
   where \f$\overline{w}_{h,\cdot}\f$ is the average weight of fish from haul \f$h\f$.
   Assume the number of fish from the haul is large so we can use the
   approximation
   \f[
   T_h\approx \frac{W_h}{E^f[w_{h,f}|\bf\theta_h]}
   \label{eq:T_h.approx}
   \f]
   where \f$\bf\theta_h\f$ is a set of parameters and random effects characterizing
   haul \f$h\f$. 
   
   Define \f$T_{a,h}\f$ to be the number of fish from haul \f$h\f$ at age \f$a\f$.
   Then 
   \f[
   T_{a,h}=T_hp_{a,h}(\bf\theta_h)
   \f]
   where \f$p_{a,h}(\bf\theta_h)\f$ is the proportion at age in haul \f$h\f$, assumed
   to be defined through the parameters and random effects characterizing
   haul \f$h\f$. Using the approximation \f$\eqref{eq:T_h.approx}\f$, we obtain
   \f[
   T_{a,h}
   \approx W_h\frac{p_{a,h}(\bf\theta_h)}{E^f[w_{h,f}|\bf\theta_h]}
   \f]
   The haul characteristics \f$\bf\theta_h\f$ are unknown, but a separate
   dataset \f$\bf D\f$ contains information of these. We consider a Bayesian
   approach where the information about \f$\bf\theta_h\f$ is given by
   the posterior \f$p(\bf\theta_h|\bf D)\f$. Assuming simulations
   \f$\bf\theta_h^s,s=1,...,S\f$ from the posterior are available, an
   estimate of \f$T_{a,h}\f$ is given by
   \f[
   \widehat T_{a,h}
   =W_h\frac{1}{S}\sum_{s=1}^S
   \frac{p_{a,h}(\bf\theta_h^s)}{E^f[w_{h,f}|\bf\theta_h^s]}
   \f]
   
   Consider now the situation where only the total weight inside a cell \f$c\f$
   is given, i.e. \f$W_c=\sum_{h\in c}W_h\f$. We could similarly as above
   define \f$p_{a,c}\f$ to be the proportion at age from cell \f$c\f$ and
   use the approximation
   \f[
   T_{a,c}
   \approx W_c\frac{p_{a,c}(\bf\theta_c)}{E^f[w_{c,f}|\bf\theta_c]}
   \f]
   Is it possible in this case only to consider simulation inside
   cells without having to simulate each haul?
   
   A problem with the use of this approximation is that
   in our modelling approach a specification of
   \f$p_{a,h}(\bf\theta_h)\f$ is given, but not 
   \f$p_{a,c}(\bf\theta_c)\f$. We do have the relationship
   \f[
   p_{a,c}=\frac{\sum_{h\in c}T_hp_{a,h}}{\sum_{h\in c}T_h}
   \f]

   where \f$n_h\f$ is the number of fish inside haul \f$h\f$. This relation 
   can not be of direct use since the \f$T_h\f$'s are unknown. 
   Note however that
   \f{eqnarray*}
   T_{a,c}
   &=&\frac{W_c}{\overline{w}_{c,\cdot}}p_{a,c}\\
   &=&\frac{W_c}{\frac{\sum_{h\in c}T_h\overline{w}_{h,\cdot}}
   {\sum_{h\in c}T_h}}
   \frac{\sum_{h\in c}T_hp_{a,h}}{\sum_{h\in c}T_h}\\
   &=&\frac{W_c}{\frac{1}{H}\sum_{h\in c}
   \frac{T_h}{\overline{T}_c}\overline{w}_{h,\cdot}}
   \frac{1}{H}\sum_{h\in c}
   \frac{T_h}{\overline{T}_c}p_{a,h}\\
   &\approx&\frac{W_c}{E^h[\frac{T_h}{\overline{T}_c}\overline{w}_{h,\cdot}]}
   E^h[\frac{T_h}{\overline{T}_c}p_{a,h}]
   \f}

   where \f$\overline{T}_c=\frac{1}{H}\sum_{h\in c}T_h\f$ and where
   \f$E^h\f$ means expectation over a random selection of hauls caught
   from cell \f$c\f$. Note that this approximation depends on a large
   number of hauls inside the cell.
   
   Assume further that \f$\frac{T_h}{\overline{T}_c}\f$ is independent from
   the other quantities involved (this assumption is questionable if
   many boats have catches close to the maximum capacity). In that case
   \f{eqnarray*}
   E^h[\frac{T_h}{\overline{T}_c}\overline{w}_{h,\cdot}]
   &=&E^h[\frac{T_h}{\overline{T}_c}]E^h[\overline{w}_{h,\cdot}]\\
   E^h[\frac{T_h}{\overline{T}_c}p_{a,h}]
   &=&E^h[\frac{T_h}{\overline{T}_c}]E^h[p_{a,h}]
   \f}
   giving
   \f[
   T_{a,c}
   \approx
   \frac{W_c}{E^h[\overline{w}_{h,\cdot}]}E^h[p_{a,h}]
   \f]
   Still, both expectation would need to be approximated over
   many Monte Carlo samples of hauls. 

   @section sec_mod The model

The model is built up in a hierarchical manner described in
Figure ??. The age \f$a\f$ is given a prior distribution
\f$p(a)\f$. The length \f$l\f$ is described through a conditional 
distribution \f$p(l|a)\f$ while finally the weight \f$w\f$ is described
through a conditional model \f$p(w|l)\f$. Note that weight given
lenght is assumed to not be influenced by age.


In the following we will define 
a cell \f$c\f$ to be a combination of year \f$y\f$, season \f$s\f$,
gear \f$g\f$ and region \f$r\f$, such that \f$c=\{y,s,g,r\}\f$. For a given \f$c\f$,
we will use \f$y(c)\f$ to denote the year corresponding to cell \f$c\f$ and so on.
 Further, let \f$h\f$ denote haul (or catch). \f$c(h)\f$ defines the cell corresponding
to haul \f$h\f$,  \f$y(h)\f$ the year corresponding to haul \f$h\f$ and so on.
The individual fish is indexed by \f$f\f$. 

Age group is indexed by \f$a\f$, where age group \f$1\f$ means the smallest age group.
\f$Y\f$ is the number of years, \f$S\f$ is the 
number of seasons, \f$G\f$ is the number of gears and \f$R\f$ is the number of regions

@subsection sec_age_mod Model for age

Define
\f{eqnarray*}
p_{h}(a)=\mbox{The proportion of fish in haul } h \mbox{ that belongs to age group } a.
\f}
We model these proportions as
\f{eqnarray*}
p_{h}(a)=\frac{\exp(\alpha_{h}^a)}{\sum_{a'}\exp(\alpha_{h}^{a'})}
\label{eq:prob.a}
\f}
where
\f{equation}
\alpha_{h}^a
=\bf X^{age}_h\bf\alpha_a+\bf Z^{age}_h\bf\zeta_a+\zeta_{h}^{age}
\f}
Here \f$\bf\alpha_a\f$ is assumed to contain fixed effects while \f$\bf\zeta_a\f$
contain random effects.

Note that there are some identifiability problems involved. In order to
solve these, we assume
\f[
\bf A\bf\alpha_a=\bf 0\\
\sum_av_a{\bf\alpha_a}=\bf 0
\f]
where \f$\bf A\f$ is a matrix defining the type of constraint. Two types of 
constraints are possible:

- Treatment constraint: In this case,
\f$\alpha^{overall,a^*}=\alpha_y^{year,a^*}=\alpha_s^{season,a^*}=\alpha_g^{gear,a^*}=0\f$
for suitable \f$a^*\f$.
Further, for identifiability of the individual \f$\alpha\f$'s, for each \f$a\f$,
we assume \f$\alpha_1^{year,a}=\alpha_1^{season,a}=\alpha_1^{gear,a}=0\f$.
We have used \f$a^*=5\f$ so far, but since the
choice of \f$a^*\f$ can effect the prior, we should look at
the sensivity in the result for different choices of \f$a^*\f$.\\
An alternative here is to make no constraints but subtract the mean after each
iteration. This should be an easy thing to modify in the program.
- Sum constraint: In this case,
\f[
\sum_a \alpha^{overall,a}=\sum_a\alpha_y^{year,a}
\f]

Again, we should look at the sensitivity of these choices compared
to using constraints on other cells.
Again we could subtract the mean after each iteration.

For the random effects, we assume
\f{eqnarray*}
[\zeta_c^{cell,a}]&=&N(0,1/\tau_{age}^{cell})\label{prior.tau.age.cell}
\f}
\f{eqnarray*}
[\zeta_{h}^{haul,a}]&=&N(0,1/\tau_{age}^{haul})\label{prior.tau.age.haul}
\f}
while \f$\{\zeta_r^{region,a},r=1,...,R \}\f$ follows an 
autoregressive CAR model,
\f{eqnarray*}
[\zeta_r^{region,a}|\zeta_{j\neq r}^{region,a}]
&=&N\left(\bar\zeta_r^{region,a},1/n_r\tau_{age}^{region}\right)
\f}
\f{eqnarray*}
\bar\zeta_r^{region,a}
&=&n_r^{-1}\sum_{j\in\delta(r)}\zeta_{j\neq r}^{region,a}
\f}
where \f$n_r\f$ is the number of neighbors of region \f$r\f$ while
\f$\delta(r)\f$ is the set of neighbors of region \f$r\f$.

For identifiability, we assume \f$\zeta_1^{region,a}=0\f$.
If \f$\tau^{region}_{age}\bf Q\f$ is the information matrix for \f$\{\zeta_r^{region,a}r=1,...,R\}\f$,
then the submatrix \f$\tau^{region}_{age}\bf Q_{-1}\f$
of \f$\tau^{region}_{age}\bf Q\f$ obtained by removing 
the first column
and the first row is the information matrix for 
\f$\{\zeta_r^{region,a}r=2,...,R\}\f$ conditional on \f$\zeta_1^{region,a}\f$.
Further, the conditional expectation is zero.


@subsection sec_mod_lga Model for length given age

Let \f$l_{h,f}\f$ be length measurement of fish \f$f\f$ from haul \f$h\f$
and \f$a_{h,f}\f$ the corresponding age. Then
\f[
[\log(l_{h,f})]=N(\beta_{0,h}+\beta_{1,h}g(a_{h,f};\bf\theta),1/\tau_{length}^{fish})
\f]
where \f$g(\cdot;\bf\theta)\f$ is a fixed function describing the relationship between age and length, while
\f$\beta_{0,h}\f$ and \f$\beta_{1,h}\f$ are haul specific intercept and slope, respectively, allowing some
variability in the relationship from one haul to another.

The function \f$g(\cdot;\bf\theta)\f$ can take many forms:

- \f$g(a)=(log(a)-log(a_{min}))\f$


The \f$\beta\f$'s are modeled as
\f{eqnarray*}
\beta_{0,h}&=&\bf X_{0,h}^{length}\bf\beta_0+\bf Z_{0,h}^{length}\bf\varepsilon_0\\
\beta_{1,h}&=&\bf X_{1,h}^{length}\bf\beta_1+\bf Z_{1,h}^{length}\bf\varepsilon_1
\f}
where \f$\bf\beta_i\f$ and \f$\bf\varepsilon_i\f$ are fixed and random effects, respectively.
\f$\bf X_{i,h}\f$ and \f$\bf Z_{i,h}\f$ are describe the covariates that influence the fixed and
random effects, respectively.
Similar constraints are needed on the fixed \f$\beta\f$-parameters:
\f{eqnarray*}
\bf B_0\bf\beta_0&=&\bf 0\\
\bf B_1\bf\beta_1&=&\bf 0
\f}
\f$\beta_1^{year}=\beta_1^{season}=\beta_1^{gear}=0\f$.
The remaining \f$\varepsilon\f$'s are random effects 
\f{eqnarray*}
[\varepsilon^{cell}_c]&=&N(0,1/\tau_{length}^{cell}),
\label{prior.tau.length.cell}
\f}
\f{eqnarray*}
[\varepsilon^{haul}_{h}]&=&N(0,1/\tau_{length}^{haul}),
\label{prior.tau.length.haul}
\f}
while \f$\{\varepsilon_r^{region,a},r=1,...,R\}\f$ follows an 
autoregressive CAR model,
\f{eqnarray*}
[\varepsilon_r^{region}|\varepsilon_{j\neq r}^{region}]
&=&N\left(\bar\varepsilon_r^{region},1/n_r\tau_{age}^{region}\right)\\
\bar\varepsilon_r^{region}
&=&n_r^{-1}\sum_{j\in\delta(r)}\varepsilon_{j\neq r}^{region}
\f}


@subsection sec_mod_wgl Model for weight given length
 
Let \f$w_{h,f}\f$ be weight measurement of fish \f$f\f$ from haul \f$h\f$. Then
\f{eqnarray*}
[\log{(w_{h,f})}]=N(\delta_{0,h} + \delta_1 \log{(l_{h,f})},1/\tau_{weight}^{fish}),
\f}
where \f$\delta_{0,h}\f$ is a cell and haul specific intercept, 
and \f$\delta_1\f$ is common to all cells and hauls. 
\f$\delta_{0,h}\f$ is given by
\f{eqnarray*}
\delta_{0,h}
= \delta^{overall} + \delta^{year}_{y(h)} + \delta^{season}_{s(h)}+
  \delta^{gear}_{g(h)}+ 
  \nu^{region}_{r(h)}+ \nu^{cell}_{c(h)}+ \nu^{haul}_{h} 
\f}
Here the \f$\delta\f$'s are fixed coefficients with constraints
\f$\delta_1^{year}=\delta_1^{season}=\delta_1^{gear}=0\f$.

The remaining \f$\nu\f$'s are
random effects
\f{eqnarray*}
[\nu^{cell}_c]&=&N(0,1/\tau_{weight}^{cell}),
\f}
\f{eqnarray*}
[\nu^{haul}_{h}]&=&N(0,1/\tau_{weight}^{haul}),
\f}
while \f$\{\nu_r^{region,a},r=1,...,R\}\f$  follows an 
autoregressive CAR model,
\f{eqnarray*}
[\nu_r^{region}|\nu_{j\neq r}^{region}]
&=&N\left(\bar\nu_r^{region},1/n_r\tau_{age}^{region}\right)\\
\bar\nu_r^{region}
&=&n_r^{-1}\sum_{j\in\delta(r)}\nu_{j\neq r}^{region}
\f}

@section sec_prior Prior distribution
The following parameters are included in the model and need 
prior distributions:
- Age parameters: \f$\alpha_y^{year,a},\alpha_s^{season,a},\alpha_g^{gear,a},\tau_{age}^{region},\tau_{age}^{cell},\tau_{age}^{haul}\f$.
- Length given age parameters: \f$\beta^{overall,a},\beta_y^{year,a},\beta_s^{season,a},\beta_g^{gear,a},\tau_{length}^{region},\tau_{length}^{cell},\tau_{length}^{haul},\tau_{length}^{fish}\f$.
- Weight given length parameters \f$\delta^{overall,a},\delta_y^{year,a},\delta_s^{season,a},\delta_g^{gear,a},\tau_{weight}^{region},\tau_{weight}^{cell},\tau_{weight}^{haul},\tau_{weight}^{fish}\f$.

We assume all \f$\alpha\f$'s, \f$\beta\f$'s and \f$\delta\f$'s are independent
and Gaussian\f$(0,K)\f$ with  \f$K\f$ large.
We further assume all \f$\tau\f$'s are Gamma\f$(B_0,B_1)\f$ with both \f$B_0\f$ and \f$B_1\f$
small.

@section sec_dat Data structure

Originally, data were assumed to be of different types:
 - Amigo data which contain complete age, length and weight data (but perhaps some missing values)
 - Reference fleet which mainly contain length-only data  but may have some additionally age measurements
   stratified by length an also some extra weight data.

   In the old version of the program, these datasets were stored differently. The Amigo data consisted
   of three vectors corresponding to age, length and weight measurements. These were ordered with respect
   to haul, with a vector nFishboat describing the number of fish in each haul. Additionally, four vectors,
   typically named f.year, f.seas, f.gear and f.area, gave the covariates corresponding to the hauls. 

   The length-only data was assumed categorised in a finite number of intervals and inside each haul, the number
   of length-measurements inside each interval was stored. The format was a Length vector storing all possible
   length-values, a Count vector giving the number of repetitions of the given length value and a Journey  vector
   giving the haul for which the observations were taken. 

   The age-stratified-by-length data was also assumed categorised in a finite number of intervals. Inside each haul,
   an ageLength vector described the length-values while ageCount was a matrix describing the numbers inside
   each age-category for that length value. A corresponding ageJourney vector gave the haul for which the
   observations were taken.

   In addition, four vectors, typically named l.year, l.seas, l.gear and l.area, gave the covariates corresponding to 
   the hauls, both for the length-only and age-stratified-by-length data. This allowed both length-only and
   age-stratified-by-length data to come from the same haul.

   For the old version, additional weight data for length-only and age-stratified-by-length data was not possible
   to use

   The new version still contain the option of these types of data. In addition, it is also opens for the possibility
   that all data is of the amigo-type, that is all data are stored in the three vectors age, length and weight.
   All values are allowed to be missing. Note though that if length is missing, only the age data is used in
   the fitting of the age model. A corresponding weight observation is currently not used. 

   \todo Utilising missing lengths better can be performed by simulating the missing length. In order to make this
   complete, the fitting of the age/lga models must be merged together with the fitting of the wgl model, which
   currently are performed through two separate calls to C-routines.
   

@section sec_ref References


\anchor Hirst_etal_2012 Hirst, D. and Storvik, G. and Rognebakke, H. and Aldrin, M. and Aanes, S. and Vølstad, J. H.
(2012)<br>
<em>A Bayesian modelling framework for the estimation of catch-at-age of commercially harvested fish species</em><br>
   Canadian Journal of Fisheries & Aquatic Sciences

\anchor Hirst_etal_2005 Hirst, D. and Storvik, G. and Aldrin, M. and Aanes, S. and Huseby, R. B. (
2005)<br>
    <em>Estimating catch-at-age by combining data from different sources</em><br>
   Canadian Journal of Fisheries & Aquatic Sciences

\anchor Hirst_etal_2004 Hirst, D. and Aanes, S. and Storvik, G. and Huseby, R. B. and Tvete, I. F.
 (2004)<br>
    <em>Estimating catch-at-age from market sampling datas using a Bayesian
        hierarchical model</em><br>
   Applied Statistics




@section sec_comp Compiling the program

The program depend on the <tt> taucs, metis, lapack and blas </tt> libraries.
   <ol>
   <li>To compile and make the caa.so file on linux, do \code R CMD SHLIB -o caa.so *.f *.c \endcode at src directory. Note that this directory should include
a <tt> Makevars </tt> file.  </li> 
   <li>To compile and make the caa.dll file on windows, we have used the MSYS evironment on windows with the
       Makefile and the dll's in the lib_win directory. With MSYS properly installed, it should be enough to
      copy all .c and .h files in the src directory together with the files in the lib_win directory to a
      common directory on a windows machine and just say \code make\endcode </li> 
  </ol>

*/
