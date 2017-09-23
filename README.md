# artool

<b><em>A</em></b>daptive <b><em>R</em></b>andomization ***tool*** for clinical trials.

The aim of the package is modelling and simulations of clinical trials with different randomiztion techniques such as

- Restricted Randomization (***RR***)
- Response-Adaptive Randomization (***RAR***)
- Covariate-Adaptive Randomization (***CA***)
- Covariate-Adjusted Response-Adaptive Randomization (***CARA***)

So far, several RR procedures have been implemented.

## Restricted Randomization procedures

Notations:
<br>
<br>
<table id="notes">
  <tr>
    <td id="center">
      <a href="https://www.codecogs.com/eqnedit.php?latex=K" target="_blank">
        <img src="https://latex.codecogs.com/gif.latex?K" title="K" />
      </a>
    </td>
    <td>Number of treatments 
      <a href="https://www.codecogs.com/eqnedit.php?latex=K&space;\geq&space;2" target="_blank">
        <img src="https://latex.codecogs.com/gif.latex?K&space;\geq&space;2" title="K \geq 2" />
      </a>
    </td>
  </tr>
  <tr>
    <td id="center"<a href="https://www.codecogs.com/eqnedit.php?latex=$w_1:w_2:\ldots:w_K$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$w_1:w_2:\ldots:w_K$" title="$w_1:w_2:\ldots:w_K$" /></a></td>
    <td>Fixed allocation ratio; $w_k$'s are positive, not necessarily equal integers with the greatest common divisor of 1</td>
  </tr>
  <tr>
    <td id="center">$\rho_k = \frac{w_k}{\sum_{k=1}^K{w_k}}$</td>
    <td>Target treatment allocation proportions $0 < \rho_k < 1$ and $\sum_{k=1}^K{\rho_k}=1$</td>
  </tr>
  <tr>
    <td id="center">$\rho=(\rho_1, \rho_2, \ldots, \rho_K)$</td>
    <td>Vector of target allocation proportions</td>
  </tr>
  <tr>
    <td id="center">$n$</td>
    <td>Total sample size for the trial</td>
  </tr>
  <tr>
    <td id="center">$N(j) = (N_1(j), N_2(j), \ldots, N_K(j))$</td>
    <td>Numbers of subjects assigned to $K$ treatments after $j$ allocations ($1 \leq j \leq n$). Note that, in general, $N_k(j)$'s are random variables with $\sum_{k=1}^K{N_k}(j)=j$</td>
  </tr>
  <tr>
    <td id="center">$P(j) = (P_{1,j}, P_{2,j}, \ldots, P_{K,j})$</td>
    <td>Vector of treatment randomization probabilities for subject $j$. Note that $0 \leq P_{k,j} \leq 1$ and $\sum_{k=1}^K{P_{k,j}}=1$ for every $j = 1, 2, \ldots, n$. Also, note that in general, $P(j)$ depends on $N(j-1)$ (in generalization of Efron's BCD) or on $\frac{N(j-1)}{j-1}$ (in generalization of Wei's UD)</td>
  </tr>
</table>
<br>
<br>


The following randomization procedures are considered:
<ol>
  <li><b>C</b>ompletely <b>R</b>andomized <b>D</b>esign (CRD): Every subject is randomizad to treatment group
  with fixed probabilities that are equal to the target allocation proportions:
  $$P_{i,j} = \rho_i, \quad i = 1, 2, \ldots, K.$$</li>

  <li><b>P</b>ermuted <b>B</b>lock <b>D</b>esign (PBD($\lambda$)): Let $W = w_1+w_2+\ldots+w_K$ and let $b =
  \lambda W$, where $\lambda$ = "number of minimal balanced sets in the block of size $b$". Let
  $k_{j-1}^* = int\left(\frac{j-1}{b}\right)$ ($int(x)$ returns the greatest integer less than or equal
  to $x$). In essence, $k_{j-1}^*$ is the number of complete blocks among the first $j-1$ assignments.
  The conditional randomization probability for the PBD($b$) is given by [Zhao and Weng (2011), page 955, equation (5)]:
  $$P_{i,j} = \frac{w_i\lambda+w_i\lambda k_{j-1}^*-N_i(j-1)}{b+b k_{j-1}^*-(j-1)},\quad i = 1, 2, \ldots, K$$</li>

  <li><b>B</b>lock <b>U</b>rn <b>D</b>esign (BUD($\lambda$)):
  This design was proposed by Zhao and Weng (2011), to provide a more random design than the PBD($b=\lambda W$). Let $N_k(j-1)$ denote the number of treatment $k$ assignments among first $j-1$ subjects, and $k_{j-1} = \text{min}_{1\leq k\leq K}\left(int\left(\frac{N_{k, j-1}}{w_k}\right)\right)$ denote the number of minimal balanced sets among the first $j-1$ assignments. Then subject $j$ is randomized to treatments with probabilities [Zhao and Weng (2011), page 955, equation (2)]:
  $$P_{i,j} = \frac{w_i\lambda+w_i k_{j-1}-N_i(j-1)}{W\lambda+W k_{j-1}-(j-1)},\quad i = 1, 2, \ldots, K$$</li>

  <li><b>M</b>ass <b>W</b>eighted <b>U</b>rn <b>D</b>esign (MWUD($\alpha$)):  This design was proposed by Zhao (2015). The parameter $\alpha > 0$ controls the maximum tolerated treatment imbalance. In his paper, Zhao (2015) considered 4 choices: $\alpha$=2; 4; 6; 8. MWUD($\alpha$) has the following formula for (conditional) treatment allocation probabilities [Zhao (2015), page 211, equation (7a)]:
  $$P_{i,j} = \frac{\text{max}\{\alpha w_i-N_i(j-1)+(j-1)w_i,0\}}{\sum_{k=1}^K\text{max}\{\alpha w_k-N_k(j-1)+(j-1)w_k,0\}},\quad i = 1, 2, \ldots, K$$
</li>

<li><b>D</b>rop-the-<b>L</b>oser Rule (DL($a$)):  The DL rule was developed by Ivanova (2003) in the context of multi-arm binary response trials with response-adaptive randomization. Here we consider its application for fixed unequal allocation. <br>
Consider an urn containing balls of $K+1$ types: types $1, 2, \ldots, K$ represent treatments and type $0$ is the immigration ball. Initially, the urn contains 1 immigration ball and $w_i$ treatment balls.  The treatment assignments are made sequentially by drawing balls at random from the urn. If a type $i$ ball is drawn ($i = 1, 2, \ldots, K$), treatment $i$ is assigned to the subject and the ball is not replaced into the urn. If type $0$ ball is drawn, it is replaced into the urn along with  $aw_1, aw_2, \ldots, aw_K$ balls of types $1, 2, \ldots, K$. The parameter $a$ is some positive integer (larger values of $a$ imply greater amount of randomness in the experiment). The DL rule is a fully randomized procedure, and it is known to be asymptotically best (Hu et al. 2006).
</li>
<br>
<li><b>D</b>oubly-Adaptive <b>B</b>iased <b>C</b>oin <b>D</b>esign (DBCD($\gamma$)): The DBCD design was developed by Hu and Zhang (2004) in the context of response-adaptive randomization. Here we consider its application for fixed unequal allocation.<br>
Initial treatment assignments ($j=1, 2, \ldots, m_0$) are made completely at random:  $P_{i,j} =\rho_i$ until each group has at least one subject (i.e. $N_k(m_0) > 0, k = 1, 2, \ldots, K$).  Subsequent treatment assignments are made as follows:
$$ P_{i,j} = \frac{\rho_i\left(\rho_i/\frac{N_i(j-1)}{j-1}\right)^\gamma}{\sum_{k=1}^K{\rho_k\left(\rho_k/\frac{N_k(j-1)}{j-1}\right)^\gamma}},\quad i = 1, 2, \ldots, K$$
where $\gamma \geq 0$ is a user-defined parameter controlling the degree of randomness ($\gamm = 0$ is most random and $\gamma\rightarrow +\infty$ is almost deterministic procedure). Importantly, treatment allocation proportions of the DBCD procedure follow an asymptotically normal distribution:

$$\sqrt{j}\left(\frac{N_i(j)}{j}-\rho_i\right) \rightarrow N\left(0, \frac{\rho_i(1-\rho_i)}{(1+2\gamma)}\right), \quad i =1, 2, \ldots, K$$
</li>

<li><b>Min</b>imum <b>Q</b>uadratic <b>D</b>istance Constrained Balance Randomization (MinQD($\eta$)):  The design was proposed by Titterington (1983), in the context of multi-arm randomized trials with covariate-adaptive randomization and balanced allocation. Here we consider an extension of this procedure to clinical trials with unequal allocation.<br>
Consider a point in the trial when $j-1$ subjects have been randomized among the $K$ treatments, and let  denote the corresponding treatment numbers ($\sum_{i=1}^K{N_i(j-1)} = j-1$). The randomization rule for the $j$th subject is as follows:
<ol type="a">
<li>For $k=1, 2, \ldots, K$, compute $B_k$, the hypothetical "lack of balance" which results from assigning the $j$th subject to treatment $K$: $B_k=\text{min}_{1\leq i\leq K}\left|\frac{N_{i,j}^{(k)}}{j}-\rho_k\right|$,  where $$N_{i,j}^{(k)} = \left\{\begin{array}{cl}N_i(j-1)+1, & i = k \\ N_i(j-1), & i\ne k\end{array}\right.$$
</li>
<li>The treatment randomization probabilities for the $j$th subject ($P_{1,j}, P_{2,j}, \ldots, P_{K,j}$) are determined as a solution to the constrained optimization problem:
$$
\begin{array}{l}
\text{minimize }\sum_{i=1}^K{\left(P_{i,j}-\rho_i\right)^2} \\
\text{subject to }\sum_{i=1}^K{B_iP_{i,j}} \leq \eta B_{(1)} + (1-\eta)\sum_{i=1}^K{B_i\rho_i} \\
\text{and }\sum_{i=1}^K{P_{i,j}}=1;\quad 0\leq P_{i,j} \leq 1, \quad i = 1,2, \ldots, K,
\end{array}
$$
where $\eta$ is the user-defined parameter (0 \leq \leq 1) that controls degree of randomness ($\eta = 0$ is most random and  $\eta = 1$ is almost deterministic procedure).
</li>
</ol>
<br>
<li><b>Max</b>imum <b>Ent</b>ropy Constraint Balance Randomization (MaxEnt($\eta$)): The design was proposed by Klotz (1978), in the context of multi-arm randomized trials with covariate-adaptive randomization and balanced allocation. Here we consider an extension of this procedure to clinical trials with unequal allocation. <br>
The MaxEnt design follows the same idea as the MinQD($\eta$) design, except for step (b), in which the constrained optimization problem deals with minimization of the Kullback-Leibler divergence between
$(P_{1,j}, P_{2,j}, \ldots, P_{K,j})$ and  $(\rho_1, \rho_2, \ldots, \rho_K)$:
$$
\begin{array}{l}
\text{minimize }\sum_{i=1}^K{P_{i,j}\text{log}\left(\frac{P_{i,j}}{\rho_i}\right)} \\
\text{subject to }\sum_{i=1}^K{B_iP_{i,j}} \leq \eta B_{(1)} + (1-\eta)\sum_{i=1}^K{B_i\rho_i} \\
\text{and }\sum_{i=1}^K{P_{i,j}}=1;\quad 0\leq P_{i,j} \leq 1, \quad i = 1,2, \ldots, K,
\end{array}
$$
where $0\leq \eta \leq 1$ controls degree of randomness of the procedure.

</li>
</ol>
<br>
<b>References</b>
<br>
<ol>
<li>Zhao W, Weng Y (2011). Block urn design—A new randomization algorithm for sequential trials with two or more treatments and balanced or unbalanced allocation. <em>Contemporary Clinical Trials</em> 32, 953-961.
</li>

<li>
Zhao W (2015). Mass weighted urn design—A new randomization algorithm for unequal allocations. <em>Contemporary Clinical Trials</em> 43, 209-216.
</li>

<li>
Ivanova A (2003). A play-the-winner-type urn design with reduced variability. <em>Metrika</em>58, 1-13.
</li>

<li>
Hu F, Zhang LX (2004). Asymptotic properties of doubly adaptive biased coin designs for multitreatment clinical trials.  <em>The Annals of Statistics</em> 32(1), 268-301.
</li>

<li>
Titterington DM (1983). On constrained balance randomization for clinical trials. <em>Biometrics</em> 39(4), 1083-1086
</li>

<li>
Klotz JH (1978). Maximum entropy constrained balance randomization in clinical trials. <em>Biometrics</em> 34(2), 283-287.
</li>
</ol>
