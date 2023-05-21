from ag import *
from qiskit.transpiler.passes import SabreSwap as router_0394   # SABRE39
from router.sabre0330_swap import SabreSwap as router_0330      # SABRE33
from router.sqgm_swap import SQGMSwap as router_sqgm            # SQGM
from router.nassc_swap import NASSCSwap as router_nassc         # NASSC


# Customise the accepted architecture graphs here
ARCHGRAPHS = {
    "ourense": ourense(),
    "tokyo": tokyo(),
    "rochester": rochester(),
    "sycamore53": sycamore53(),
    "sycamore54": sycamore54()
}

# Customise the accepted routers here
ROUTERS = {
    "sabre0330": router_0330,
    "sabre0394": router_0394,
    "sqgm": router_sqgm,
    "nassc": router_nassc
}