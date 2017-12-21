# Math 315 Cheat Sheet, By Andre Kaba

## ***First order***

$$y'=ky \implies \mathbf{y=ce^{kx}}$$

---

### **Graphical**

$$y'=F(x,y)$$
1) Draw (x,y) axis
1) Draw at some points a small line with slope $F(x,y)$$
1) Draw Integral curves tangent to the lines (each one is a solution)
1) Integral curves don't intersect cuz $F(x,y)$ only have one value

![Example](integralCurves.png)

### **Numerical**

$N$: number of steps

$T$: we are approximating on $[0,T]$

$h=\frac{T}{N}$: approximation step

Forward euler:

$$\hat{X}(t+h)=X(t)+hf(X(t))$$
Backward euler:

$$\hat{X}(t+h)=X(t)+hf(X(t+h))$$
Trapazoid:

$$\hat{X}(t+h)=X(t)+\frac{h}{2}(f(X(t))+f(X(t+h)))$$
Improved euler:

$$\hat{X_{*}}=X(t)+hf(X(t))$$
$$\hat{X}(t+h)=X(t)+\frac{h}{2}\left(f(X(t))+f(\hat{X_{*}})\right)$$
Also,

LTE= $|X(t_1)-\hat{X}(t_1)|=h^2\ or\ h^3$

GTE= $N*LTE \sim T h^{\alpha-1}$

#### **Matlab code:**

##### Improved euler

```matlab
clear all;
T=1;
N=100;
h=T/N;

t=0;
x=1;

for ts=1:N
    t=t+h;
    xstart=x+h*x;
    x=x+0.5*h*(x+xstar);
    xstore(t)=x;
    tstore(t)=t;
end
xexact=exp(tstore);

plot(tstore,xexact,'ro');
hold on;
plot(tstore,xstore,'bx');

legend('exact','improved');
xlabel('time');
ylabel('x(t)');
```

### **Seperable**

$$N(y)\frac{dy}{dx}=M(x)$$
Solution: integrate both sides.

### **Linear**

Standard form:

$$y'+p(t)y=q(t)$$
Homogeneous if $q(t)=0$ and seperable. Otherwise we get the solution by integrating factors:

$$u(t)=e^{\int p(t) dt}$$
$$y(t)=\frac{\int u(t)q(t)dt +c}{u(t)}$$

#### Tempreture Example

$$\dot{T}=k(T_{ext}-T)$$
Seperable if $T_{ext}$ is constant

### **Substitution $y'=F(\frac{y}{x})$**

Becomes seperable by substitution:

$$Z(x)=\frac{y}{x} \implies y=xz,\ y'=z+xz'$$
$$\implies \frac{dz}{F(z)-z}=\frac{dx}{z}$$

### **Exact Equations**

We have the form:

$$M(x,y)+N(x,y)\frac{dy}{dx}=0$$
We need to find a function

$$\Psi(x,y): \Psi_x=M(x,y),\Psi_y=N(x,y)$$

Then $\Psi(x,y)=c$

1) Test $M_y=N_x$
1) Integrate either $\int M dx=...+h(y)$ or $\int N dy=...+h(x)$ to get $\Psi(x,y)$
1) Use the other integral to find $h$ using a seperable ODE (usually)
 
$$N=...+h'(y)\ or\ M=...+h'(x)$$
1) Set $\Psi(x,y)=....+h(x)\ or\ h(y) +k$ to get implcit solution
1) Simplify further if needed

#### Converting to exact

We multiply by integrating factor $U$ to get:

$$UM+UN\frac{dy}{dx}=0$$
We want:

$$(UM)_y=(UN)_x \implies MU_y+M_y-NU_x-N_x=0$$

We have three cases:

1) assume $U$ is a function of just $x$, $\implies U_Y=0 \implies \frac{dU}{dx}=\frac{M_y-N_x}{N}U$

1) assume $U$ is a function of just $y$, $\implies U_X=0 \implies \frac{dU}{dx}=\frac{N_X-M_Y}{M}U$

1) we can try $U(x,y)=x^a y^b \implies U_x=ax^{a-1}y^b,U_y=x^aby^{b-1}$ and sub in the original equation to see what might happen.

Probably You will get a seperable equation in the end to solve for $U$.

### **Bernoulli**

We have the form:

$$y'+p(x)y=q(x)y^n$$

1) divde by $y^n$
1) Sub $z=y^{1-n}$
1) $z'=(1-n)y^{-n}y'$
1) plug in the original equation to get a linear one

$$\frac{1}{1-n}z'+p(x)z=q(x)$$
which is linear!

---

\pagebreak

## ***Second order***

### **Linear with constant coefficients**

We have the form

$$ay''+by'+cy=q(t)$$
When $q(t)=0$ we call it homogeneuos

#### Homogenous:

---

##### Superposition principle

If $y_1(t)$ and $y_2(t)$ are two solutions to the linear homogenueous equation then so:

$$y(t)=c_1y_1(t)+c_2y_2(t)$$
This work with non constant coeffs too!

All solutions to the Homogenous equations are of the form $e^{rt}$ where we get $r$ from the charactristic equation:

$$ar^2+br+c=0$$

And we have 3 cases:
1) **Real roots $r_1\neq r_2$ :**

    Solution is $y(t)=c_1e^{r_1t}+c_2e^{r_2t}$

1) **Complex roots $r_{_{1,2}}=a\pm ib$ :**

    Solution is $y(t)=c_{_1}e^{at}cos(bt)+c_{_2}e^{at}sin(bt)$

    In addition, if $Z$ is a complex-valued solution, then $Re(Z)$ and $Im(Z)$ are both lineary independent solutions too!
    **Note:** proof for euler formula and relationship between $cos,sin,e$ all start with the number $Z=cos(t)+isin(t)$

1) **Double roots $r_1= r_2=r$ :**

    Solution is $y(t)=c_{_1}e^{rt}+c_{_2}te^{rt}$
---

#### Homogenous but without constant coeffs (Reduction of order):

We have:

$$p(t)y''+q(t)y'+r(t)y=0$$
$$y_1(t)=f(t)$$

We can generate another solution $y_2(t)$ and get the general solution by superposition principal as follows:

1) Take $y_2(t)=v(t)y_1(t)$
1) Sub $y_2,y_2',y_2''$ in the original equation you should lose the $v$ term to get $f(t)v''+g(t)v'=0$
1) substitute $w=v'$ then solve the linear equation for $w$ and get $v=\int w dt$
1) take any constants in $v$ to simplify it's most general form
1) then compute $y_2(t)$ and use superposition to get a general solution

---

#### **Non-Homogenous special cases:**

##### 1) **ERF (when RHS is cos,sin,e...):**

We have:

$$ay''+by'+cy=F_{ext}(t)$$


We can do something when $F_{ext}(t)$ is $sin,cos,e,e*cos,e*sin...$

**Note:** sometimes $F_{ext}$ can be refered to as input signal

**Strategy:**
1) Find particular solution to the equation $y_p$
1) Find the homogeneuos solution $y_h(t)$ (assume $F_{ext}=0$)
1) General solution is $y(t)=y_p(t)+y_h(t)$

In general, input signal is an exponential, $ay''+by'+cy=Ae^{rt}$ and we use the guess:$y_p(t)=Be^{rt}$ where $B=\frac{A}{P(r)}; P(r)\neq0$ as a particular solution for the ODE. Here $P(s)=as^2+bs+c$ is the characteristic poly of the ODE.
So the particular solution can be obtained fro ERF forumla:

$$y_p=\frac{A}{P(r)}e^{rt}$$

**Notice 1** that when $P(r)=0$ we have what's called resonance.

**Notice 2:**
If the RHS is $cos/sin$ then we can complexify the solution to get $e^{it}$ in the RHS then user ERF

We can generlize the ERF if $P(r)=0$:

* If $P(r)=0,P'(r)=0...,P^{(m-1)}(r)=0$ but $P^{(m)}(r)\neq0$ then we can use the formula to compute a particular solution $y_p(t)$

$$y_p(t)=\frac{At^me^{rt}}{P^{(m)}(r)}$$

---

##### 2) **Undeterminal coefficients (when input signal is a poly):**

We have:

$$P(D)x=b_kt^k+...+b_1t+b_0$$
$$P(0)\neq0$$

The solution is a poly of degree $\leq k$, we assume it's $X(t)=a_kt^k+...+a_1t+a_0$ then derivate and sub in the ODE and match with RHS to compute the $a_i$

**Note:** If $P(0)=0$ we can use the substituation $u=y'$ to reduce order and solve as linear first order ODE.

---

##### 3) **Variation of parameters(when $A$ is not a constant in $...=Ae^{rt}$**

We have

$$P(D)x=A(t)e^{rt}$$

Solution: assume we have a solution $y_p(t)=U(t)e^{rt}$ then derivate and sub in the original ODE  to get an ODE solvable by Undeterminal coefficients.

---
\pagebreak

## ***Fourier Transform***

We have a periodic function $f(t)=f(t+2L)$ which can be harmonic too i.e periodic for $2L,4L,6L,...$. So any window $[-L,L]$ determine the function.

* Linear combination of periodic functions is also periodic: $af(t)+bg(t)$ is periodic for example.

We define fourier series as the linear combination of all harmonics:

$$f(t)=\frac{a_0}{2}+a_1cos(t)+a_2cos(2t)+a_3cos(3t)+...$$
$$+b_1sin(t)+b_2sin(2t)+b_3sin(3t)+...$$
Where:

$$L=period[-\pi,\pi]/2$$

$$a_0=\frac{1}{L} \int_{-L}^{L} f(t)dt$$

$$a_n =\frac{1}{L} \int_{-L}^{L} f(t)cos(nt)dt$$

$$b_n= \frac{1}{L} \int_{-L}^{L} f(t)sin(nt)dt$$

So we can use this to represent the RHS of some ODE with $sin/cos$ and we know how to solve it in these cases.

* Note: to show a function is continous at $a_0$ we must have:

    $$\lim{x \rightarrow a_0^{-}} f(x)=\lim{x \rightarrow a_0^{+}} f(x)=f(a_0)$$
### Useful formulas

$$\int (Ax^2+Bx+c)cos(n\pi x)dx=\frac{\pi n(2Ax+B)cos(n\pi x)+sin(n\pi x)\left[A(\pi^2 n^2 x^2-2)+\pi^2n^2(Bx+c)\right]}{\pi^3 n^3}+c_2$$

$$\int (Ax^2+Bx+c)sin(n\pi x)dx=\frac{\pi n(2Ax+B)sin(n\pi x)-cos(n\pi x)\left[A(\pi^2 n^2 x^2-2)+\pi^2n^2(Bx+c)\right]}{\pi^3 n^3}+c_2$$

$$sin(2t)=2sin(t)cos(t)$$

$$sin(A)cos(B)=\frac{sin(A+B)-sin(A-B)}{2}$$

$$sin^2(t)=\frac{1-cos(2t)}{2}$$

### Odd/even function

$$f \text{ is even } \implies f(-x)=f(x)$$
$$f \text{ is odd } \implies f(-x)=-f(x)$$

## **Laplace transform**

### Solving IVP using Laplace

We solve ODE using laplace transform by moving to Laplace space and solving there algebraically then doing inverse transform to get $y(t)$.
Use this fact for Laplace of derivatives:

$$\mathbb{L}[f^{(n)}]=s^nF(s)-s^{n-1}f(0)-s^{n-2}f'(0)-...-sf^{n-2}(0) -f^{n-1}(0)$$

Extra notes to solve problems:

* initial conditions should always start from 0, if not (it start from $k$) transform them to 0 by using $\nu=t-k$ and $u(\nu)=y(\nu+k)$ and solve for $\nu$ then for $y$

#### ![Laplace table](Laplacetable.png)

#### Rules for partial fractions

When you see a demoninator of the following forms split them as follows:

$$ ax+b  \implies \frac{A}{ax+b}$$
$$ (ax+b)^k  \implies \frac{A_1}{ax+b}+\frac{A_2}{(ax+b)^2}+...+\frac{A_k}{(ax+b)^k}$$

$$ ax^2+bx+c  \implies \frac{Ax+B}{ax^2+bx+c}$$
$$ (ax^2+bx+c)^k  \implies \frac{A_1x+B_1}{ax^2+bx+c}+\frac{A_2x+B_2}{(ax^2+bx+c)^2}+...+\frac{A_kx+B_k}{(ax^2+bx+c)^k}$$

### Step function $u_c(t)=u(t-c)$

#### ![step function <>](stepfunction.png)

It have these properties:

$$\mathbb{L}[u(t-c)f(t-c)]=e^{-cs}F(s)$$

taking $f(t-c)=1$ we get:

$$\mathbb{L}[u(t-c)]=\frac{e^{-cs}}{s}$$

### Dirac delta function

used when a force is highly affecting the system in small amount of time. $\delta(t-a)$ has the following properties:

$$\delta(t-a)=0, t \neq a$$

$$\int_{a-\epsilon}^{a+\epsilon} \delta(t-a) dt=1,\epsilon >0 $$

$$\int_{a-\epsilon}^{a+\epsilon} f(t)\delta(t-a) dt=f(t),\epsilon >0$$

$$\mathbb{L}[\delta(t-a)]\int_{0}^{\infty} e^{-at}\delta(t-a) dt=e^{-as},a >0$$

### Convolutions

Some times we might face a situation were we want the $\mathbb{L}^{-1}[F(s)G(s)]$. Hence we need a convolution, i.e convulation in real space is multiplication in laplace space!
We need to remember this: (* is convolution)

$$(f*g)(t)=(g*f)(t)=\int_{0}^{t} f(t-\tau)g(\tau)d\tau=\int_{0}^{t} f(\tau)g(t-\tau)d\tau$$

$$\mathbb{L}[(f*g)(t)]=F(s)G(s)$$
$$\mathbb{L}^{-1}[F(s)G(s)]=(f*g)(t)$$

## **System of equations**

### Converting $2^{nd}$ order equation to system of two equations:

1) let $x_1(t)=y(t)$, $x_2(t)=y'(t)$
1) derivate $x_1'(t)=y'(t)=x_2(t)$ and $x_2'(t)=y''(t)=$ get something from original equation in term of $x_1(t),x_2(t)$
1) Compute initial conditions $x_1(t_0)=y(t_0)$, $x_2(t_1)=y'(t_1)$

### Solving $\vec{X}'=A \vec{X}$

1) Eigen values

    $$det(A-\lambda I)=0$$
1) Eigen vectors

    $$(A-\lambda_k I)\vec{v_{k}}=\vec{0}$$
    for $k=1,2,3...$
1) Solution is:
    1) Real unique Eigen values: $X_g(t)=c_1\vec{v_1}e^{\lambda_1 t}+c_2\vec{v_2}e^{\lambda_2 t}$
    1) Complex Eigen values: if they have imaginary parts and real parts we get a spiral outwards if real part is positive and inwards. Otherwise, if only imaginary we get an ellipse centered at the origin. The general solution would be

    $$x(t)=c_1u(t)+c_2v(t)$$

    where we get $u(t),v(t)$ from the first solution $x_1(t)=u(t)+iv(t)$
    1) Double root (degenerate node): we need to compute the first eigen vector $u(t)$ then find another vector $(A-\lambda I)\rho=u$ then the general solution is

    $$x_1(t)=c_1e^{\lambda t}u+c_2(te^{\lambda t}u+e^{\lambda t} \rho)$$

### Solving non linear shit

We have

$$x'=f(x,y)$$
$$y'=g(x,y)$$
where $f,g$ might be non linear!

Solution steps:
1) find a fixed point (or more) $(x^*,y^*)$ that makes $f,g$ equal to 0 at it.
1) Do change of variables to it $u=x-x^*,v=y-y^*$ (no need appearently)
1) compute the new jacobian $(A')$

$$
J(x,y)=
\left(\begin{array}{cc}
f_x & f_y\\
g_x & g_y
\end{array}\right)
$$

1) solve for the new jacobian $u=Jv$

## **Power series**

### First type

We want to solve:

$P(x)y''+Q(x)y'+R(x)y=0$

Solution:
1) let $y(x)=\sum_{n=1}^{\infty} c_nx^n$
1) find $y',y''$ and sub in the original equation
1) Solve for $c_n$

### Second type (frobenious series crap)

We want to solve:

$y''+P(x)y'+Q(x)y=0$ with a regular singular point $x=a$ that makes $P \text{ or } Q \rightarrow \infty$ as $x \rightarrow a$

Solution:
1) let $p(x)=xP(x), q(x)=x^2Q(x)$ and rewrite the equation to $x^2y''+p(x)y'+q(x)y=0$
1) Check if both $p(x),q(x)$ are analytical on $x=a$ i.e they posses all derivatives and can be expanded into a powe series (they should be. Otherwise, $x=a$ is irregular singular point)
1) write $y(x)=\sum_{n=0}^{\infty} a_nx^{n+r}$ and sub in the original equation
1) Solve the indecial equation $cr(r-1)+p_0r+q_0=0$ (we get it from the coefficients of $a_0$) and pick the larger root $r$
1) Sub in the original equation and solve for $a_n$
1) solution would be $y(x)=x^r\sum_{n=0}^{\infty} a_nx^{n}$

## Digression on linear independence

We have:

$$y''+p(x)y'+q(x)y=0$$
We can generlize to any order IVP:

$$y^{(n)}+p_{n-1}(x)y^{(n-1}+p_{n-2}(x)y^{(n-2)}+...+p_0(x)=0$$

Solution:
1) compute the wronskion (can be used to check independence if $\neq 0$), $y_i$ is a solution

$$W=det\left(
\left(\begin{array}{cc}
y_1 & y_2\\
y_1' & y_2'
\end{array}\right)
    \right)$$

1) solve the following new problem $W'+p(x)w=0$ using abel's identity

$$w(t)=w_0e^{-\int_{t_0}^{t}p(t)dt},w_0=w(t_0)$$

## Eular equations (one possible view point)

we have an equation of this form:

$$x^2y''+axy'+by=0$$

solution:
1) assume we have solution $y(x)=x^m$ (if we assume $y(x)=e^u$ we get same shit as chapter two)
1) subtitute in the original equation and compute the roots of the charactristic poly
1)  * Two real roots: $y(x)=c_1x^{m_1}+c_2x^{m_2}$
    * double root: $y(x)=c_1x^{m}\ln(x)+c_2x^{m}$
    * complex root: $y(x)=c_1x^\alpha \cos(\beta\ln{x})+c_2x^\alpha \sin(\beta\ln{x})$
