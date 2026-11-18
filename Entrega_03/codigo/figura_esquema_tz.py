"""Esquema que explica la corrección del tiempo a vertical, tz = t z / R."""
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]   # carpeta Entrega_03
import matplotlib.pyplot as plt, numpy as np
fig,ax=plt.subplots(figsize=(10,7.5))
ax.fill_between([-1.5,6.5],[0,0],[-8,-8],color="#f1e6c8")
ax.plot([-1.5,6.5],[0,0],"k",lw=2)
ax.add_patch(plt.Rectangle((-0.15,-7.5),0.3,7.5,color="#d0d0d0",ec="k"))
ax.text(0,0.35,"Pozo",ha="center",fontsize=12)
ax.plot(3,0,"v",ms=16,color="#c0392b"); ax.text(3,0.35,"Fuente (golpe)",ha="center",fontsize=12,color="#c0392b")
z=5
ax.plot(0,-z,"o",ms=12,color="#1f4e79"); ax.text(-0.35,-z,"Geófono\n(canal 5)",ha="right",va="center",fontsize=11,color="#1f4e79")
ax.annotate("",xy=(0.15,-z+0.1),xytext=(3,0),arrowprops=dict(arrowstyle="->",lw=3,color="#c0392b"))
ax.text(1.85,-2.7,"R = √(z² + d²)\nla onda viaja en diagonal\ntarda t (lo que se lee)",fontsize=12,color="#c0392b")
ax.annotate("",xy=(0.55,-z),xytext=(0.55,0),arrowprops=dict(arrowstyle="<->",lw=2,color="#1f4e79"))
ax.text(0.7,-4.3,"z (vertical)",fontsize=12,color="#1f4e79",rotation=90)
ax.annotate("",xy=(3,-0.6),xytext=(0,-0.6),arrowprops=dict(arrowstyle="<->",lw=1.5))
ax.text(1.5,-0.95,"d = 3 m",ha="center",fontsize=12)
ax.text(3.6,-5.0,"Idea: si la onda hubiera bajado derecho\nrecorrería z en vez de R.\nSe escala el tiempo con la misma proporción:",fontsize=12)
ax.text(3.6,-6.3,"tz = t · z / R",fontsize=20,weight="bold",color="#1f4e79")
ax.text(3.6,-7.4,"Canal 5: z = 5 m, R = √(25+9) = 5,83 m,\nt = 41 ms  →  tz = 41·5/5,83 = 35,2 ms",fontsize=11.5)
ax.set_xlim(-1.5,9.5); ax.set_ylim(-8,1); ax.axis("off")
fig.tight_layout(); fig.savefig(ROOT / "figuras" / "esquema_tz.png", dpi=200)
