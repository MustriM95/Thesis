# Microbial Consumer Resource Model 

$ \Large \frac{1}{N_{i}} \frac{dN_{i}}{dt} = \sum \limits _{\alpha = 1} ^{M} \Delta w_{i \alpha} C_{i\alpha} R_{\alpha} - m_{i} $

$ \Large \frac{dR_{\beta}}{dt} = S_{\beta} - \sum \limits _{i = 1} ^{S} C_{i\beta} R_{\beta} N_{i} + \sum \limits _{\alpha = 1} ^{M} \sum \limits _{i=1} ^{S}
D_{\beta \alpha}^{i}C_{i \alpha} N_{i} R_{\alpha}  $

| Symbol | Definition | Units |
| :- | :- | :-: |
| $N_{i}$ | Biomass of $i^{th}$ consumer | mass |
| $ R_{\alpha} $ | Concentration of the $\alpha^{th}$ resource | mass/volume |
| $ \Delta w_{i \alpha} $ | Biomass produced by species $i$ given uptake of a volume of resource $\alpha$ after metabolite leakage | mass/mass/volume |
| $ C_{i \alpha} $ | Uptake rate of species $i$ per unit concentration of resource $\alpha$ | 1/mass \*time |
| $ m_{i} $ |  Maitenance biomass utilization of $i^{th}$ species | 1/time |
| $ S_{\beta} $ | Supply of resource $\beta$ | mass/volume \*time |
| $ D_{\alpha \beta }^{i} $ | Stoichiometric matrix, number of molecules of resource $\beta$ secreted by species $i$ through the consumption of resource $\alpha$ | unitless |

Additionally $\Delta w_{i\alpha} = w_{\alpha} - \sum \limits _{\beta = 1} ^{M} D_{\beta \alpha}^{i}w_{\beta} $.
Where $w_{\alpha}$ is the ATP yield of resource $\alpha$ and $\sum \limits _{\beta = 1} ^{M} D_{\beta \alpha}^{i}w_{\beta}$ is the total ATP yield of all metabolic byproducts.

# Constraints

## Energy conservation
To guarantee that no energy/mass is created through the leaking of metabolites we impose a simple constraint on $\Delta w_{i \alpha} $. Stated explicitely, the maximum ATP yeild of resource $\alpha$ must exceed the total ATP yield of all metabolic byproducts. 
<br>
$ w_{\alpha} > \sum \limits _{\beta = 1} ^{M} D_{\beta \alpha}^{i} w_{\beta} $

## Metabolic tradeoffs
From a biological standpoint, microbes can only allocate a finite amount of cellular resources to nutrient uptake and are thus limited in the quantity and diversity of nutrients they are able to utilize. This can be integrated into the model by constraining the structure of the nutrient uptake matrix $C_{i \alpha}$. While there are a variety of ways to accomplish this, the original constraint takes into account non-randomly distributed consumers sharing similar consumption coefficients between taxonomically related groups. 


To begin with, family-level consumption parameters are deterimined by two parameters $\theta_{\alpha, f}$, the concentration parameter of resource $\alpha$ belonging to family $f$, and $\Omega_{f}$, the magnitude of all concentration parameters such that $ \Omega_{f} = \sum \limits _{\alpha} \theta_{\alpha, f} $.


```julia

```
