# Regenera las figuras de los apuntes:
#   uv run --with matplotlib,numpy figs/generar_figuras.py
# (o bien: python3 figs/generar_figuras.py, con matplotlib y numpy instalados)
#
# Cada figura reproduce exactamente los datos de las tablas de los capítulos;
# los asserts al final comprueban que la numerología no se ha desviado.

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGS = os.path.dirname(os.path.abspath(__file__))


def fig_euler_comparacion():
    # 01-PVI: y' = -2y, y(0)=1, h=0.1 — exacta vs Euler adelante/atrás
    h = 0.1
    x = np.arange(0, 1 + h / 2, h)
    xf = np.linspace(0, 1, 400)
    y_ad = (1 - 2 * h) ** (x / h)
    y_at = (1 + 2 * h) ** (-x / h)
    fig, ax = plt.subplots(figsize=(6, 3.6))
    ax.plot(xf, np.exp(-2 * xf), "k-", lw=2, label=r"Exacta $e^{-2x}$")
    ax.plot(x, y_ad, "o--", color="tab:blue", label="Euler adelante")
    ax.plot(x, y_at, "s--", color="tab:red", label="Euler atrás")
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
    ax.set_title(r"$y'=-2y$, $h=0.1$: la exacta queda encuadrada")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{FIGS}/euler_comparacion.png", dpi=150)
    plt.close(fig)
    return abs(y_ad[5] - 0.3277) < 1e-3 and abs(y_at[5] - 0.4019) < 1e-3


def fig_euler_pasos():
    # 01-PVI: ejemplo paso a paso — Euler adelante sobre y'=-2y con h=0.25
    # (paso grueso a proposito para que se vea la escalera de tangentes):
    # cada punto nuevo cae sobre la recta de pendiente -2*y_n y el rayo se
    # prolonga un tramo mas alla para mostrar donde seguiria esa direccion
    h = 0.25
    n = np.arange(5)
    x = n * h
    y = (1 - 2 * h) ** n                      # 1, 0.5, 0.25, 0.125, 0.0625
    xf = np.linspace(0, 1.05, 400)
    fig, ax = plt.subplots(figsize=(6, 3.6))
    ax.plot(xf, np.exp(-2 * xf), "k-", lw=2, label=r"Exacta $e^{-2x}$")
    ax.plot(x, y, "o", color="tab:blue", ms=7, label="Euler adelante")
    for i in range(4):
        xs = np.linspace(x[i] - 0.03, x[i] + 1.35 * h, 2)
        ax.plot(xs, y[i] - 2 * y[i] * (xs - x[i]), ":", color="tab:red",
                lw=1.5, alpha=0.8,
                label=r"tangente $f=-2y_n$" if i == 0 else None)
   
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
    ax.set_title(r"Euler adelante paso a paso: $y'=-2y$, $h=0{,}25$")
    ax.legend(loc="upper right"); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{FIGS}/euler_pasos.png", dpi=150)
    plt.close(fig)
    return abs(y[4] - 0.0625) < 1e-9 and abs(y[1] - 0.5) < 1e-9


def fig_rk_convergencia():
    # 02-RungeKutta: error global en x=1 vs h, pendientes 1, 2, 4
    def err(method, N):
        h = 1.0 / N
        y = 1.0
        f = lambda y: -2 * y
        for _ in range(N):
            if method == "euler":
                y += h * f(y)
            elif method == "rk2":
                k1 = f(y); k2 = f(y + h / 2 * k1)
                y += h * k2
            else:
                k1 = f(y); k2 = f(y + h / 2 * k1)
                k3 = f(y + h / 2 * k2); k4 = f(y + h * k3)
                y += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return abs(y - np.exp(-2))
    Ns = np.array([4, 8, 16, 32, 64, 128])
    hs = 1.0 / Ns
    exacta = np.exp(-2)
    fig, ax = plt.subplots(figsize=(6, 3.6))
    for nombre, metodo, c in [("Euler ($h$)", "euler", "tab:blue"),
                              ("RK2 ($h^2$)", "rk2", "tab:orange"),
                              ("RK4 ($h^4$)", "rk4", "tab:green")]:
        e = [err(metodo, int(n)) for n in Ns]
        ax.loglog(hs, e, "o-", color=c, label=nombre)
    for p, c in [(1, "tab:blue"), (2, "tab:orange"), (4, "tab:green")]:
        ax.loglog(hs, hs ** p * 0.15, ":", color=c, alpha=0.6)
    ax.set_xlabel("$h$"); ax.set_ylabel("error global en $x=1$")
    ax.set_title(r"Convergencia sobre $y'=-2y$: pendientes $h$, $h^2$, $h^4$")
    ax.legend(); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(f"{FIGS}/rk_convergencia.png", dpi=150)
    plt.close(fig)
    # un paso RK4 con h=0.1 debe dar 0.8187333 como en el texto
    h, y, f = 0.1, 1.0, lambda y: -2 * y
    k1 = f(y); k2 = f(y + h / 2 * k1); k3 = f(y + h / 2 * k2); k4 = f(y + h * k3)
    return abs(y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4) - 0.8187333) < 1e-6


def fig_rk4_esquema():
    # 02-RungeKutta: un paso RK4 real sobre y' = 2 - y, y0 = 0.5, h = 1.
    # Las cuatro pendientes de prueba k_i (violeta) se evaluan en los puntos
    # huecos que alcanzan los rayos punteados; el avance real es una unica
    # recta gruesa roja con la pendiente ponderada phi — un paso tipo Euler
    # con esa pendiente — que cae casi sobre la exacta 2 - 1.5 e^{-x}
    f = lambda x, y: 2.0 - y
    x0, y0, h = 0.0, 0.5, 1.0
    k1 = f(x0, y0)
    p2 = (x0 + h / 2, y0 + h / 2 * k1);  k2 = f(*p2)
    p3 = (x0 + h / 2, y0 + h / 2 * k2);  k3 = f(*p3)
    p4 = (x0 + h, y0 + h * k3);          k4 = f(*p4)
    phi = (k1 + 2 * k2 + 2 * k3 + k4) / 6
    y1 = y0 + h * phi                     # k = 1.5,.75,1.125,.375; y1 = 1.4375
    yex = lambda x: 2 - 1.5 * np.exp(-x)  # exacta de y' = 2 - y, y(0)=0.5

    fig, ax = plt.subplots(figsize=(6.4, 4))
    xf = np.linspace(-0.02, 1.12, 300)
    ax.plot(xf, yex(xf), "k-", lw=2, label=r"exacta $2-1.5e^{-x}$")
    # rayos predictores: del punto inicial a cada punto de prueba
    for px, py in (p2, p3, p4):
        ax.plot([x0, px], [y0, py], ":", color="tab:pink", lw=1.2)
    # pendientes de prueba: segmento corto de pendiente k_i en su punto
    for (px, py), kk, nom, off in zip(
            [(x0, y0), p2, p3, p4], [k1, k2, k3, k4],
            ["$k_1$", "$k_2$", "$k_3$", "$k_4$"],
            [(-40, -5), (-35, 0), (-35, -30), (-35, 5)]):
        d = 0.10
        ax.plot([px - d, px + d], [py - d * kk, py + d * kk],
                color="tab:purple", lw=2.2)
        ax.annotate(nom, (px + d, py + d * kk), textcoords="offset points",
                    xytext=off, color="tab:purple", fontsize=13)
    # puntos: inicio lleno, prueba huecos, llegada roja sobre la exacta
    ax.plot(x0, y0, "o", color="tab:blue", ms=7)
    ax.annotate(r"$y_n$", (x0, y0), textcoords="offset points",
                xytext=(6, -13), color="tab:red", fontsize=14)
    for px, py in (p2, p3, p4):
        ax.plot(px, py, "o", ms=6, mfc="white", mec="tab:purple")
    ax.plot(x0 + h, y1, "o", color="tab:red", ms=7)
    ax.annotate(r"$y_{n+1}$", (x0 + h, y1), textcoords="offset points",
                xytext=(-15, -15), color="tab:red", fontsize=14)
    # el avance real: una sola recta con la pendiente ponderada phi
    ax.plot([x0, x0 + h], [y0, y1], color="tab:red", lw=2.8,
            label=r"paso real: pendiente $\varphi$")
    ax.annotate(r"$\varphi=\frac{k_1+2k_2+2k_3+k_4}{6}$", (0.70, 1.00),
                color="tab:red", fontsize=14)
    # guias verticales en el punto medio y el final del paso
    for xv, lab in [(x0, r"$x_n$"), (x0 + h / 2, r"$x_n+h/2$"), (x0 + h, r"$x_n+h$")]:
        ax.axvline(xv, color="darkgray", lw=0.9, ls=":")
        ax.annotate(lab, (xv, 0.18), ha="center", fontsize=14,
                    color="gray")
    ax.set_xlim(-0.08, 1.16); ax.set_ylim(0.1, 1.95)
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
    ax.set_title("Un paso de RK4: cuatro pendientes y una ponderada",
                 fontsize=11)
    ax.legend(loc="upper left", fontsize=10); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{FIGS}/rk4_esquema.png", dpi=150)
    plt.close(fig)
    return (abs(y1 - 1.4375) < 1e-9
            and abs(yex(1) - y1) < 0.02)   # y1=1.4375 vs exacta ~1.4482


def fig_rk2_esquema():
    # 02-RungeKutta: un paso RK2 real sobre y' = -2y, y0 = 1, h = 0.5
    # (h grueso para que se vea la geometria; el ejemplo del texto usa 0.1).
    # k1 pendiente en el inicio; el rayo punteado es el medio paso Euler hasta
    # el punto medio hueco; k2 es la pendiente ahi; el avance real es una
    # unica recta gruesa roja con pendiente k2 desde (x_n, y_n). El rombo gris
    # marca donde hubiera caido Euler puro con el mismo h (cruza el eje).
    f = lambda x, y: -2.0 * y
    x0, y0, h = 0.0, 1.0, 0.5
    k1 = f(x0, y0)
    pm = (x0 + h / 2, y0 + h / 2 * k1);  k2 = f(*pm)
    y1 = y0 + h * k2                      # k1=-2, k2=-1, y1=0.5
    y_euler = y0 + h * k1                 # Euler puro con h: 0
    yex = np.exp

    fig, ax = plt.subplots(figsize=(6.4, 4))
    xf = np.linspace(-0.02, 0.62, 300)
    ax.plot(xf, yex(-2 * xf), "k-", lw=2.6, label=r"exacta $e^{-2x}$")
    # predictor: medio paso Euler con k1, y prolongado hasta x_n+h (Euler)
    ax.plot([x0, x0 + h], [y0, y_euler], ":", color="tab:pink", lw=1.6)
    ax.plot(x0 + h, y_euler, "D", ms=6, mfc="white", mec="gray")
    ax.annotate("Euler", (x0 + h, y_euler), textcoords="offset points",
                xytext=(5, -4), color="gray", fontsize=12)
    # pendientes de prueba k1 (inicio) y k2 (punto medio)
    for (px, py), kk, nom, off in zip(
            [(x0, y0), pm], [k1, k2], ["$k_1$", "$k_2$"],
            [(-32, -10), (-35, -10)]):
        d = 0.06
        ax.plot([px - d, px + d], [py - d * kk, py + d * kk],
                color="tab:purple", lw=3)
        ax.annotate(nom, (px + d, py + d * kk), textcoords="offset points",
                    xytext=off, color="tab:purple", fontsize=13)
    # puntos: inicio lleno, medio hueco, llegada roja
    ax.plot(x0, y0, "o", color="tab:blue", ms=7)
    ax.annotate(r"$y_n$", (x0, y0), textcoords="offset points",
                xytext=(6, 0), color="tab:red", fontsize=13)
    ax.plot(*pm, "o", ms=6, mfc="white", mec="tab:purple")
    ax.plot(x0 + h, y1, "o", color="tab:red", ms=7)
    # avance real: una sola recta con la pendiente del punto medio
    ax.plot([x0, x0 + h], [y0, y1], color="tab:red", lw=3.6,
            label=r"paso real: pendiente $k_2$")
    ax.annotate(r"$y_{n+1}=y_n+h\,k_2$", (x0 + h, y1),
                textcoords="offset points", xytext=(5, 5),
                ha="left", color="tab:red", fontsize=13)
    for xv, lab in [(x0, r"$x_n$"), (x0 + h / 2, r"$x_n+h/2$"), (x0 + h, r"$x_n+h$")]:
        ax.axvline(xv, color="gray", lw=0.9, ls=":")
        ax.annotate(lab, (xv, -0.14), ha="center", fontsize=13,
                    color="dimgray")
    ax.set_xlim(-0.05, 0.66); ax.set_ylim(-0.18, 1.25)
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
    ax.set_title("Un paso de RK2: la pendiente del punto medio", fontsize=13)
    ax.legend(loc="upper right", fontsize=11); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{FIGS}/rk2_esquema.png", dpi=150)
    plt.close(fig)
    return (abs(y1 - 0.5) < 1e-9
            and abs(y_euler) < 1e-9)       # Euler puro cae en 0 con h=0.5


def fig_cc_ejemplo():
    # 04-CondicionesContorno: misma EDO y'' = -1 (barra con fuente de calor
    # uniforme), tres juegos de condiciones. Dirichlet-Dirichlet da el arco
    # simetrico; la mixta (borde derecho aislado, y'=0) la asimetrica con
    # tangente horizontal; Neumann puro no tiene solucion porque el calor
    # generado no puede evacuar por ningun extremo.
    x = np.linspace(0, 1, 200)
    yD = 0.5 * x * (1 - x)                  # y(0)=y(1)=0
    yM = x - 0.5 * x**2                     # y(0)=0, y'(1)=0
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.plot(x, yD, color="tab:blue", lw=2.5,
            label=r"Dirichlet: $y(0)=y(1)=0$")
    ax.plot(x, yM, color="tab:orange", lw=2.5,
            label=r"Mixta: $y(0)=0$, $y'(1)=0$")
    # marcas de borde: punto = valor fijo, tangente horizontal = aislado
    ax.plot(0, 0, "o", color="k", ms=6)
    ax.plot(1, 0, "o", color="tab:blue", ms=6)
    ax.annotate("$y=0$", (1, 0), textcoords="offset points",
                xytext=(-6, -16), ha="right", fontsize=10)
    ax.plot(1, 0.5, "o", color="tab:orange", ms=6)
    ax.plot([0.8, 1.1], [0.5, 0.5], ":", color="tab:orange", lw=2)
    ax.annotate("$y'=0$: aislado", (1, 0.5), textcoords="offset points",
                xytext=(-8, 10), ha="right", color="tab:orange", fontsize=10)
    ax.set_xlim(-0.05, 1.15); ax.set_ylim(-0.08, 0.62)
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
    ax.set_title(r"Misma EDO $y''=-1$, distintas condiciones de contorno",
                 fontsize=11)
    ax.legend(loc="upper left", fontsize=9); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{FIGS}/cc_ejemplo.png", dpi=150)
    plt.close(fig)
    return (abs(yD[-1]) < 1e-12                # y(1)=0 en la Dirichlet
            and abs((1 - 1.0)) < 1e-12)        # y'(1)=0 en la mixta


def fig_euler_estabilidad():
    # 03-Estabilidad: y' = -10y, Euler adelante con lambda*h = 1.5, 2 y 2.5;
    # dos paneles porque la divergente (hasta ~25) aplastaria a las estables
    lam = 10.0
    xf = np.linspace(0, 2, 400)
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.4), sharex=True)
    paneles = [(axes[0], [(0.15, "tab:blue"), (0.2, "tab:green")], 1.6,
                r"Estable ($\lambda h=1.5$) y marginal ($\lambda h=2$)"),
               (axes[1], [(0.25, "tab:red")], 30,
                r"Inestable ($\lambda h=2.5$): diverge")]
    for ax, casos, ylim, titulo in paneles:
        ax.plot(xf, np.exp(-lam * xf), "k-", lw=2,
                label=r"Exacta $e^{-10x}$")
        for h, c in casos:
            n = int(round(2 / h))
            x = np.arange(n + 1) * h
            ax.plot(x, (1 - lam * h) ** np.arange(n + 1), "o--", color=c,
                    ms=4, label=rf"Euler adelante, $\lambda h = {lam*h:.1f}$")
        ax.set_xlabel("$x$"); ax.set_ylabel("$y$"); ax.set_ylim(-ylim, ylim)
        ax.set_title(titulo); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{FIGS}/euler_estabilidad.png", dpi=150)
    plt.close(fig)
    # los tres regimenes: |1-lam*h| = 0.5, 1.0, 1.5
    g = np.abs(1 - lam * np.array([0.15, 0.2, 0.25]))
    return g[0] < 1 and g[1] == 1 and g[2] > 1


def fig_euler_rigida():
    # 03-Estabilidad: reacciones consecutivas A -(k1=100)-> B -(k2=1)-> C.
    # El intermediario y2 contiene los dos modos en una sola funcion; si el
    # rapido diverge, el acoplamiento k1*h*y1 arrastra la medida.
    # (a) [0, 0.35]: adelante con h=0.021 (k1*h = 2.1, 5% sobre el limite)
    # explota; atras con h=0.1 sigue. (b) [0, 1]: adelante a su paso legal
    # h=0.019 (52 pasos, con oscilacion moribunda) vs atras h=0.1 (10 pasos)
    k1, k2 = 100.0, 1.0
    Bex = lambda x: (k1 / (k1 - k2)) * (np.exp(-k2 * x) - np.exp(-k1 * x))

    def integrar(h, n, atras):
        y1 = np.empty(n + 1); y2 = np.empty(n + 1)
        y1[0], y2[0] = 1.0, 0.0
        for i in range(n):
            if atras:   # (I - hA) y_{n+1} = y_n, triangular: sustitucion
                y1[i + 1] = y1[i] / (1 + k1 * h)
                y2[i + 1] = (y2[i] + k1 * h * y1[i + 1]) / (1 + k2 * h)
            else:
                y1[i + 1] = (1 - k1 * h) * y1[i]
                y2[i + 1] = y2[i] + h * (k1 * y1[i] - k2 * y2[i])
        return y2

    # paneles apilados: el HTML muestra img a max-width 36em (~576px); una
    # figura 1425x510 se reducia a ~40%, ilegible
    fig, axes = plt.subplots(2, 1, figsize=(6.4, 6.4))
    # (a) 5% sobre el limite: la oscilacion de y1 contamina y2
    ax = axes[0]
    h = 0.021
    n = int(round(0.35 / h))                   # 17 pasos, x_max = 0.357
    x = np.arange(n + 1) * h
    xf = np.linspace(0, x[-1], 300)
    ax.plot(xf, Bex(xf), "k-", lw=2.4, label=r"Exacta $y_2(x)$")
    ax.plot(x, integrar(h, n, False), "o--", color="tab:blue", ms=5,
            label=rf"Euler adelante, $h={h}$")
    xb = np.arange(4) * 0.1
    ax.plot(xb, integrar(0.1, 3, True), "s--", color="tab:red", ms=6,
            label=r"Euler atrás, $h=0.1$")
    ax.set_xlabel("$x$"); ax.set_ylabel("$y_2$"); ax.set_ylim(-6.5, 6.5)
    ax.set_title(r"5% sobre el límite: el modo muerto contamina $y_2$",
                 fontsize=11)
    ax.legend(fontsize=9.5); ax.grid(alpha=0.3)
    # (b) paso legal de adelante vs atras
    ax = axes[1]
    xf = np.linspace(0, 1, 300)
    ax.plot(xf, Bex(xf), "k-", lw=2.4, label=r"Exacta $y_2(x)$")
    h = 0.019
    n = int(1 / h)                             # 52 pasos, x_max = 0.988
    x = np.arange(n + 1) * h
    ax.plot(x, integrar(h, n, False), "o-", color="tab:blue", ms=3, lw=1.2,
            label=rf"Euler adelante, $h={h}$")
    xb = np.arange(11) * 0.1
    ax.plot(xb, integrar(0.1, 10, True), "s--", color="tab:red", ms=6,
            label=r"Euler atrás, $h=0.1$")
    ax.set_xlabel("$x$"); ax.set_ylabel("$y_2$"); ax.set_ylim(-0.6, 2)
    ax.set_title(r"$y_2$ en $[0,1]$: adelante necesita 52 pasos, atrás 10",
                 fontsize=11)
    ax.legend(fontsize=9.5); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{FIGS}/euler_rigida.png", dpi=150)
    plt.close(fig)
    return (abs(integrar(0.021, 17, False)[-1]) > 3         # explota
            and abs(integrar(0.1, 10, True)[-1] - Bex(1)) / Bex(1) < 0.1)


def fig_disparo():
    # 05-PVF: y''=y, trayectorias s*sinh(x) hacia y(1)=1
    x = np.linspace(0, 1, 200)
    s_star = 1 / np.sinh(1)
    fig, ax = plt.subplots(figsize=(6, 3.6))
    for s, c, lab in [(0.5, "tab:blue", "$s_1=0.5$ (corto)"),
                      (1.3, "tab:orange", "$s_2=1.3$ (se pasa)"),
                      (s_star, "tab:green", rf"$s^*={s_star:.4f}$ (da en $\beta$)")]:
        ax.plot(x, s * np.sinh(x), "--", color=c, label=lab)
    ax.axhline(1, color="k", lw=1, ls=":")
    ax.plot([1], [1], "ko")
    ax.annotate(r"$\beta = 1$", (1, 1), textcoords="offset points",
                xytext=(-30, 8))
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
    ax.set_title(r"Método del disparo sobre $y''=y$, $y(0)=0$, $y(1)=1$")
    ax.legend(loc="upper left"); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{FIGS}/disparo.png", dpi=150)
    plt.close(fig)
    return abs(s_star * np.sinh(1) - 1) < 1e-12


def disparo_psi(lam, h=0.005):
    # 06-Autovalores: integra psi''=-lam*psi, psi(0)=0, psi'(0)=1 con RK4
    # sobre [0,1]; devuelve (x, psi). Nada usa la solucion cerrada.
    n = int(round(1 / h))
    x = np.linspace(0, 1, n + 1)
    u = np.array([0.0, 1.0])              # u = (psi, psi')
    f = lambda v: np.array([v[1], -lam * v[0]])
    psi = np.empty(n + 1)
    psi[0] = 0.0
    for i in range(n):
        k1 = f(u); k2 = f(u + h / 2 * k1)
        k3 = f(u + h / 2 * k2); k4 = f(u + h * k3)
        u = u + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        psi[i + 1] = u[0]
    return x, psi


def fig_autovalores_F():
    # 06-Autovalores: barrido REAL — cada punto F(lam)=psi(1) sale de integrar
    # el PVI con RK4; las raices se refinan por biseccion sobre el integrador
    lams = np.arange(1.0, 101.0, 1.0)
    F = np.array([disparo_psi(l)[1][-1] for l in lams])
    raices = []
    for i in range(len(lams) - 1):
        if F[i] * F[i + 1] >= 0:
            continue
        lo, hi = lams[i], lams[i + 1]
        flo = F[i]
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            fm = disparo_psi(mid)[1][-1]
            if flo * fm <= 0:
                hi = mid
            else:
                lo, flo = mid, fm
        raices.append(0.5 * (lo + hi))
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    ax.plot(lams, F, "o", ms=3, color="tab:blue", alpha=0.85)
    ax.axhline(0, color="k", lw=0.8)
    ax.axvspan(9, 16, color="tab:red", alpha=0.12,
               label="cambio de signo $\\lambda\\in[9,16]$")
    for n, r in enumerate(raices, 1):
        ax.plot([r], [0], "o", color="tab:red")
        ax.annotate(rf"$\lambda_{n}$", (r, 0), textcoords="offset points",
                    xytext=(0, -14), ha="center")
    ax.set_xlabel(r"$\lambda$"); ax.set_ylabel(r"$F(\lambda)$")
    ax.set_title(r"Barrido del residuo $F(\lambda)=\psi(1;\lambda)$: "
                 "cada punto integra el PVI")
    ax.legend(loc="upper right"); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{FIGS}/autovalores_F.png", dpi=150)
    plt.close(fig)
    return len(raices) == 3 and abs(raices[0] - np.pi ** 2) < 1e-6


def fig_autovalores_disparo():
    # 06-Autovalores: trayectorias integradas con RK4 (el disparo real);
    # solo lam = pi^2 cierra en psi(1)=0
    fig, ax = plt.subplots(figsize=(6, 3.6))
    for lam, c, lab in [(5.0, "tab:blue", r"$\lambda=5$: llega por encima"),
                        (np.pi ** 2, "tab:green", r"$\lambda=\pi^2$: da en 0"),
                        (14.0, "tab:orange", r"$\lambda=14$: pasa de largo")]:
        x, psi = disparo_psi(lam)
        ax.plot(x, psi, "--", color=c, label=lab)
        ax.plot([1], [psi[-1]], "o", color=c)
    ax.axhline(0, color="k", lw=0.8)
    ax.axvline(1, color="gray", lw=0.8)
    ax.annotate(r"$\psi(1)=0$", (1, 0), textcoords="offset points",
                xytext=(-52, -16))
    ax.set_xlabel("$x$"); ax.set_ylabel(r"$\psi$")
    ax.set_title(r"Disparo en $\lambda$: solo $\lambda=\pi^2$ cierra en $\psi(1)=0$")
    ax.legend(loc="upper left", fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{FIGS}/autovalores_disparo.png", dpi=150)
    plt.close(fig)
    x, psi = disparo_psi(np.pi ** 2)
    return abs(psi[-1]) < 1e-6


def fig_laplace_placa():
    # 08-Elipticas: Laplace en [0,1]^2, borde superior 100, resto 0 (Gauss-Seidel)
    N = 64
    u = np.zeros((N + 1, N + 1))
    u[-1, :] = 100.0
    for _ in range(6000):
        u[1:-1, 1:-1] = 0.25 * (u[2:, 1:-1] + u[:-2, 1:-1]
                                + u[1:-1, 2:] + u[1:-1, :-2])
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    im = ax.imshow(u, origin="lower", extent=[0, 1, 0, 1], cmap="hot")
    fig.colorbar(im, ax=ax, label="temperatura")
    for i, (xj, yi) in enumerate([(1 / 3, 1 / 3), (2 / 3, 1 / 3),
                                 (1 / 3, 2 / 3), (2 / 3, 2 / 3)], 1):
        ax.plot(xj, yi, "o", color="cyan", ms=7, mfc="none")
        ax.annotate(f"$u_{i}$", (xj, yi), textcoords="offset points",
                    xytext=(6, 6), color="cyan")
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
    ax.set_title("Ecuación de Laplace: borde superior a 100, resto a 0")
    fig.tight_layout(); fig.savefig(f"{FIGS}/laplace_placa.png", dpi=150)
    plt.close(fig)
    k = N // 3
    # malla fina: los nodos del ejemplo grueso (12.5, 37.5) dan ~11.6, ~36.5
    return 10 < u[k, k] < 13 and 34 < u[2 * k, k] < 39


def _idx(c, d):
    return c if d == 0 else f"{c}{d:+d}"


def _fig_malla(con_stencil, nombre, paso_j="h"):
    # grilla 3x3 generica con etiquetas u_{i,j} y pasos h, paso_j;
    # con_stencil resalta la cruz que conecta u_{i,j} con sus 4 vecinos
    fig, ax = plt.subplots(figsize=(2.9, 2.7))
    g = [-1, 0, 1]
    for v in g:
        ax.axhline(v, color="gray", lw=0.8)
        ax.axvline(v, color="gray", lw=0.8)
    if con_stencil:
        # cruz del stencil
        for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
            ax.plot([0, dx], [0, dy], color="tab:green", lw=3, alpha=0.6,
                    zorder=2)
    for x in g:
        for y in g:
            if con_stencil:
                color = ("tab:red" if (x, y) == (0, 0)
                         else "tab:blue" if x == 0 or y == 0 else "silver")
                tinte = "k" if x == 0 or y == 0 else "gray"
            else:
                color, tinte = "tab:blue", "k"
            ax.plot(x, y, "o", ms=9, zorder=3, color=color)
            ax.annotate(rf"$u_{{{_idx('i',x)},{_idx('j',y)}}}$", (x, y),
                        textcoords="offset points", xytext=(8, 6),
                        fontsize=7.5, color=tinte)
    # separacion h en ambas direcciones
    ax.annotate("", (0, -1.28), (-1, -1.28),
                arrowprops=dict(arrowstyle="<->", color="dimgray", lw=1))
    ax.annotate("$h$", (-0.5, -1.28), ha="center", va="top",
                textcoords="offset points", xytext=(0, -4), fontsize=10)
    ax.annotate("", (-1.28, 0), (-1.28, -1),
                arrowprops=dict(arrowstyle="<->", color="dimgray", lw=1))
    ax.annotate(f"${paso_j}$", (-1.28, -0.5), ha="right", va="center",
                textcoords="offset points", xytext=(-4, 0), fontsize=10)
    ax.annotate("$x$", (1.52, -1.45), ha="right", fontsize=11)
    ax.annotate("$y$", (-1.45, 1.5), ha="left", fontsize=11)
    ax.set_xlim(-1.7, 1.7); ax.set_ylim(-1.7, 1.7)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"Malla uniforme: $x_i=ih$, $y_j=j{paso_j}$", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIGS}/{nombre}.png", dpi=150)
    plt.close(fig)
    return True


def fig_malla_notacion():
    # 07-EDP: grilla uniforme e indices u_{i,j} — convencion comun a 08 y 09
    # (paso k en la segunda direccion; en las elipticas k = h)
    return _fig_malla(False, "malla_notacion", "k")


def fig_stencil_laplace():
    # 08-Elipticas: la misma grilla; el stencil conecta u_{i,j} con sus
    # 4 vecinos — las esquinas no participan
    return _fig_malla(True, "stencil_laplace")


def _fig_stencil_xt(nombre, cn):
    # rejilla x-t 3x2 (pasos h, k); cn=False dibuja la T de FTCS (los tres
    # nodos de la capa j alimentan a u_{i,j+1}), cn=True el bloque de
    # Crank-Nicolson (las tres incognitas de j+1 se resuelven acopladas)
    fig, ax = plt.subplots(figsize=(3.4, 2.5))
    gx, gt = [-1, 0, 1], [0, 1]
    for v in gx:
        ax.axvline(v, color="gray", lw=0.8)
    for v in gt:
        ax.axhline(v, color="gray", lw=0.8)
    # acoples del stencil
    if cn:
        for dx in gx:
            ax.plot([dx, dx], [0, 1], color="tab:green", lw=3, alpha=0.6,
                    zorder=2)
        for t in gt:
            ax.plot([-1, 1], [t, t], color="tab:green", lw=3, alpha=0.6,
                    zorder=2)
    else:
        for dx in gx:
            ax.plot([dx, 0], [0, 1], color="tab:green", lw=3, alpha=0.6,
                    zorder=2)
    for x in gx:
        for t in gt:
            en = cn or t == 0 or x == 0
            color = ("tab:red" if t == 1 and (cn or x == 0)
                     else "tab:blue" if en else "silver")
            ax.plot(x, t, "o", ms=9, zorder=3, color=color)
            ax.annotate(rf"$u_{{{_idx('i',x)},{_idx('j',t)}}}$", (x, t),
                        textcoords="offset points", xytext=(8, 6),
                        fontsize=7.5, color="k" if en else "gray")
    # separacion h (espacio) y k (tiempo)
    ax.annotate("", (0, -0.28), (-1, -0.28),
                arrowprops=dict(arrowstyle="<->", color="dimgray", lw=1))
    ax.annotate("$h$", (-0.5, -0.28), ha="center", va="top",
                textcoords="offset points", xytext=(0, -4), fontsize=10)
    ax.annotate("", (-1.35, 1), (-1.35, 0),
                arrowprops=dict(arrowstyle="<->", color="dimgray", lw=1))
    ax.annotate("$k$", (-1.35, 0.5), ha="right", va="center",
                textcoords="offset points", xytext=(-4, 0), fontsize=10)
    ax.annotate("$x$", (1.62, -0.52), ha="right", fontsize=11)
    ax.annotate("$t$", (-1.55, 1.42), ha="left", fontsize=11)
    ax.set_xlim(-1.75, 1.75); ax.set_ylim(-0.7, 1.65)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"Stencil {'Crank–Nicolson' if cn else 'FTCS'}: "
                 f"$x_i=ih$, $t_j=jk$", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIGS}/{nombre}.png", dpi=150)
    plt.close(fig)
    return True


def fig_stencil_ftcs():
    # 09-Parabolicas: T espacio-tiempo del esquema explicito
    return _fig_stencil_xt("stencil_ftcs", cn=False)


def fig_stencil_cn():
    # 09-Parabolicas: bloque 3x2 del esquema implicito
    return _fig_stencil_xt("stencil_cn", cn=True)


def fig_malla_laplace():
    # 08-Elipticas: malla 4x4 con contorno conocido e incognitas u1..u4
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    g = [0, 1 / 3, 2 / 3, 1]
    for v in g:
        ax.axhline(v, color="gray", lw=0.8)
        ax.axvline(v, color="gray", lw=0.8)
    for x in g:
        for y in g:
            borde = x in (0, 1) or y in (0, 1)
            if borde:
                val = 100 if y == 1 else 0
                ax.plot(x, y, "s", ms=9,
                        color="tab:red" if val else "tab:blue")
                ax.annotate(f"{val}", (x, y), textcoords="offset points",
                            xytext=(9, 4), fontsize=8, color="dimgray")
    nodos = {"u_1": (1 / 3, 1 / 3), "u_2": (2 / 3, 1 / 3),
             "u_3": (1 / 3, 2 / 3), "u_4": (2 / 3, 2 / 3)}
    for nombre, (x, y) in nodos.items():
        ax.plot(x, y, "o", ms=13, mfc="white", color="k", mew=1.5)
        ax.annotate(rf"${nombre}$", (x, y), ha="center", va="center",
                    fontsize=9)
    # stencil de 5 puntos sobre u1
    cx, cy = nodos["u_1"]
    for dx, dy in [(1 / 3, 0), (0, 1 / 3), (0, -1 / 3)]:
        ax.annotate("", (cx + dx, cy + dy), (cx, cy),
                    arrowprops=dict(arrowstyle="-", color="tab:green",
                                    lw=2, alpha=0.7))
    ax.annotate("", (cx - 1 / 3, cy), (cx, cy),
                arrowprops=dict(arrowstyle="-", color="tab:green",
                                lw=2, alpha=0.7))
    ax.set_xticks(g); ax.set_yticks(g)
    ax.set_xticklabels(["0", "1/3", "2/3", "1"])
    ax.set_yticklabels(["0", "1/3", "2/3", "1"])
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
    ax.set_title("Malla $h=1/3$: contorno conocido, 4 incógnitas interiores")
    ax.set_xlim(-0.18, 1.18); ax.set_ylim(-0.18, 1.18)
    ax.set_aspect("equal")
    fig.tight_layout(); fig.savefig(f"{FIGS}/malla_laplace.png", dpi=150)
    plt.close(fig)
    return True


def _difusion_superficie(nombre, Un, k, etiqueta):
    # superficie exacta u(x,t) translucida + solucion numerica sobre la
    # malla real (h=0.25) como malla facetada
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    nlev = Un.shape[0]
    x = np.linspace(0, 1, 60)
    t = np.linspace(0, (nlev - 1) * k, 40)
    X, T = np.meshgrid(x, t)
    U = np.exp(-np.pi ** 2 * T) * np.sin(np.pi * X)
    fig = plt.figure(figsize=(6.8, 4.6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, T, U, cmap="hot", alpha=0.45, linewidth=0)
    xm = np.array([0, 0.25, 0.5, 0.75, 1.0])
    tn = np.arange(nlev) * k
    for j, tj in enumerate(tn):
        ax.plot(xm, np.full(5, tj), Un[j], "-", color="gray", lw=0.8)
        ax.plot(xm, np.full(5, tj), Un[j], "o", color="gray", ms=4,
                label=etiqueta if j == 0 else None)
    for i in range(5):
        ax.plot(np.full(nlev, xm[i]), tn, Un[:, i], "-", color="gray",
                lw=0.8)
    # condicion inicial: curva t=0; contorno: bordes x=0 y x=1
    ax.plot(x, np.zeros_like(x), np.sin(np.pi * x), "r-", lw=3,
            label="condición inicial")
    for xv, lab in [(0, "contorno"), (1, None)]:
        ax.plot([xv] * len(t), t, np.zeros(len(t)), "b-", lw=3, label=lab)
    ax.set_xlabel("$x$"); ax.set_ylabel("$t$"); ax.set_zlabel("$u$")
    ax.view_init(elev=22, azim=-62)
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title(rf"$u_t=u_{{xx}}$: solución exacta (superficie) "
                 rf"y {etiqueta} (malla)", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{FIGS}/{nombre}.png", dpi=150)
    plt.close(fig)


def fig_difusion_superficie():
    # FTCS r=1/2: u_{i,j+1} = 0.5*(u_{i-1,j} + u_{i+1,j})
    k, nlev = 0.03125, 6          # r = k/h^2 = 1/2; niveles 0 .. 0.15625
    Un = np.zeros((nlev, 5))
    Un[0] = np.sin(np.pi * np.array([0, 0.25, 0.5, 0.75, 1.0]))
    for j in range(1, nlev):
        Un[j, 1:4] = 0.5 * (Un[j - 1, 0:3] + Un[j - 1, 2:5])
    _difusion_superficie("difusion_superficie", Un, k, "FTCS $r=1/2$")
    return abs(Un[2, 2] - 0.5) < 1e-12 and abs(Un[4, 1] - 0.1768) < 1e-3


def fig_difusion_superficie_cn():
    # misma vista con Crank-Nicolson r=1 (k=0.0625, 4 niveles: 0 .. 0.1875)
    k, nlev = 0.0625, 4
    A = np.array([[2, -0.5, 0], [-0.5, 2, -0.5], [0, -0.5, 2]])
    Un = np.zeros((nlev, 5))
    Un[0] = np.sin(np.pi * np.array([0, 0.25, 0.5, 0.75, 1.0]))
    for j in range(1, nlev):
        b = 0.5 * (Un[j - 1, 0:3] + Un[j - 1, 2:5])
        Un[j, 1:4] = np.linalg.solve(A, b)
    _difusion_superficie("difusion_superficie_cn", Un, k, "CN $r=1$")
    return abs(Un[3, 2] - 0.1636) < 5e-3


def fig_difusion_ftcs_cn():
    # 09-Parabolicas: u_t=u_xx, h=0.25 — FTCS r=1 (inestable) | FTCS r=1/2 | CN r=1
    h, D = 0.25, 1.0
    x = np.array([0.25, 0.5, 0.75])
    xf = np.linspace(0, 1, 200)

    def ftcs_paso(u, r):
        uL = np.concatenate(([0.0], u[:-1]))   # borde u = 0
        uR = np.concatenate((u[1:], [0.0]))
        return u + r * (uL + uR - 2 * u)

    def perfil(ax, u, t, lab):
        linea, = ax.plot([0, *x, 1], [0, *u, 0], "o-", ms=4, label=lab)
        ax.plot(xf, np.exp(-np.pi ** 2 * t) * np.sin(np.pi * xf), ":",
                color=linea.get_color(), lw=2, alpha=0.9)

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.4), sharey=True)

    # (a) FTCS r=1, condicion inicial con perturbacion de alta frecuencia;
    # la exacta (punteada) incluye el modo 3 que decae como exp(-9 pi^2 t)
    u = np.sin(np.pi * x) + 0.05 * np.sin(3 * np.pi * x)
    ax = axes[0]
    for n in range(4):
        linea, = ax.plot([0, *x, 1], [0, *u, 0], "o-", ms=4, label=f"$j={n}$")
        t = n * 0.0625
        uex = (np.exp(-np.pi ** 2 * t) * np.sin(np.pi * xf)
               + 0.05 * np.exp(-9 * np.pi ** 2 * t) * np.sin(3 * np.pi * xf))
        ax.plot(xf, uex, ":", color=linea.get_color(), lw=2, alpha=0.9)
        if n < 3:
            u = ftcs_paso(u, 1.0)
    ax.plot([], [], ":", color="gray", lw=2, label="exacta")
    ax.set_title(r"FTCS, $r=1$ (inestable)")
    ax.set_xlabel("$x$"); ax.set_ylabel("$u$")
    ax.set_ylim(-0.5, 1)
    ax.legend(fontsize=6); ax.grid(alpha=0.3)
    ok_a = np.allclose(u, [-0.4472, 0.7746, -0.4472], atol=1e-3)

    # (b) FTCS r=1/2, condicion inicial suave (k = 0.03125)
    u = np.sin(np.pi * x)
    ax = axes[1]
    perfil(ax, u, 0.0, "$t=0$")
    for n in range(1, 5):
        u = ftcs_paso(u, 0.5)
        if n in (2, 4):
            perfil(ax, u, n * 0.03125, rf"$t={n*0.03125:g}$")
    ax.plot([], [], ":", color="gray", lw=2, label="exacta")
    ax.set_title(r"FTCS, $r=1/2$ (estable)")
    ax.set_xlabel("$x$")
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    ok_b = np.allclose(u, [0.1768, 0.25, 0.1768], atol=1e-3)

    # (c) Crank-Nicolson r=1 (k = 0.0625)
    r = 1.0
    A = np.array([[1 + r, -r / 2, 0], [-r / 2, 1 + r, -r / 2], [0, -r / 2, 1 + r]])
    u = np.sin(np.pi * x)
    ax = axes[2]
    perfil(ax, u, 0.0, "$t=0$")
    for n in range(1, 4):
        uL = np.concatenate(([0.0], u[:-1]))
        uR = np.concatenate((u[1:], [0.0]))
        b = r / 2 * (uL + uR) + (1 - r) * u
        u = np.linalg.solve(A, b)
        perfil(ax, u, n * 0.0625, rf"$t={n*0.0625:g}$")
    ax.plot([], [], ":", color="gray", lw=2, label="exacta")
    ax.set_title(r"Crank–Nicolson, $r=1$")
    ax.set_xlabel("$x$")
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    ok_c = abs(u[1] - 0.1636) < 5e-3  # CN u_2 en t=0.1875 segun la tabla

    fig.tight_layout(); fig.savefig(f"{FIGS}/difusion_ftcs_cn.png", dpi=150)
    plt.close(fig)
    return ok_a and ok_b and ok_c


if __name__ == "__main__":
    checks = [
        ("euler_comparacion", fig_euler_comparacion),
        ("euler_pasos", fig_euler_pasos),
        ("rk_convergencia", fig_rk_convergencia),
        ("rk4_esquema", fig_rk4_esquema),
        ("rk2_esquema", fig_rk2_esquema),
        ("cc_ejemplo", fig_cc_ejemplo),
        ("euler_estabilidad", fig_euler_estabilidad),
        ("euler_rigida", fig_euler_rigida),
        ("disparo", fig_disparo),
        ("autovalores_F", fig_autovalores_F),
        ("autovalores_disparo", fig_autovalores_disparo),
        ("laplace_placa", fig_laplace_placa),
        ("malla_notacion", fig_malla_notacion),
        ("stencil_laplace", fig_stencil_laplace),
        ("malla_laplace", fig_malla_laplace),
        ("stencil_ftcs", fig_stencil_ftcs),
        ("stencil_cn", fig_stencil_cn),
        ("difusion_superficie", fig_difusion_superficie),
        ("difusion_superficie_cn", fig_difusion_superficie_cn),
        ("difusion_ftcs_cn", fig_difusion_ftcs_cn),
    ]
    for nombre, fn in checks:
        ok = fn()
        print(f"{'OK ' if ok else 'FALLA'} {nombre}")
        assert ok, nombre
    print("Figuras generadas en", FIGS)
