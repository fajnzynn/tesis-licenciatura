"""
=============================================================================
PASO 1 — OBTENER CONTACTOS ENTRE RESIDUOS
=============================================================================

OBJETIVO:
    Obtener una lista de pares de residuos que están en contacto en cada
    estado conformacional (abierto y cerrado), junto con la fracción del
    tiempo que ese contacto está formado (TotalFrac).

HAY DOS RUTAS POSIBLES. Elegí la que corresponda a tus datos:

  RUTA A — Desde PDB cristalográfico (sin trayectoria)
  ─────────────────────────────────────────────────────
  Usala cuando:
    • Solo tenés el archivo .pdb descargado del Protein Data Bank
    • Querés testear el pipeline con una proteína nueva rápidamente
    • No tenés simulaciones all-atom previas del sistema

  Qué hace: calcula distancias euclidianas entre todos los pares de Cα
  en la estructura estática y define como contacto toda distancia < 8 Å.
  Como la estructura es estática, TotalFrac = 1 para todos los contactos
  (o se filtra directamente por distancia).

  RUTA B — Desde trayectoria all-atom (cpptraj)
  ───────────────────────────────────────────────
  Usala cuando:
    • Tenés simulaciones MD all-atom (.nc + .parm7) del sistema
    • Querés contactos promediados sobre el ensemble conformacional
    • Necesitás TotalFrac realista (fracción de frames con contacto)

  Qué hace: cpptraj analiza la trayectoria y calcula para cada par (i,j)
  qué fracción del tiempo los Cα están a menos de 8 Å.

OUTPUT COMÚN DE AMBAS RUTAS:
    res_contactsAbierto.dat    →  contactos del estado abierto
    res_contactsCerrado.dat    →  contactos del estado cerrado

    Formato:
        #Res1   #Res2   TotalFrac   Contacts
        1       5       0.87        1
        1       6       0.45        1
        ...

    #Res1 y #Res2 son numeración base 1 (igual que PDB).

DEPENDENCIAS:
    RUTA A: biopython, pandas, numpy
    RUTA B: cpptraj (AMBER), luego este script para leer los resultados

=============================================================================
"""

import pandas as pd
import numpy as np
from itertools import combinations

# =============================================================================
# PARÁMETROS — editá estos valores para tu sistema
# =============================================================================

# Archivos PDB de cada estado (RUTA A)
PDB_ABIERTO  = "1z98_assembly.pdb"    # estructura del estado abierto (ej: Protein Data Bank)
PDB_CERRADO  = "2b5f.pdb"    # estructura del estado cerrado

# Archivos de trayectoria (RUTA B)
PARM_ABIERTO  = "alfacarbonsAbierto.parm7"
TRAJ_ABIERTO  = "alfacarbonsAbierto.nc"
PARM_CERRADO  = "alfacarbonsCerrado.parm7"
TRAJ_CERRADO  = "alfacarbonsCerrado.nc"

# Parámetros del análisis de contactos
DISTANCIA_CONTACTO_A  = 8.0   # Å — umbral para definir contacto (RUTA A)
DISTANCIA_CONTACTO_B  = 8     # Å — igual pero para cpptraj (RUTA B)
OFFSET_SECUENCIA      = 4     # ignorar pares con |i-j| < este valor
                              # (vecinos en secuencia ya cubiertos por enlace/ángulo/diedro)

# Outputs
OUT_ABIERTO  = "res_contactsAbierto.dat"
OUT_CERRADO  = "res_contactsCerrado.dat"


# =============================================================================
# RUTA A — Contactos desde PDB cristalográfico (usando Biopython)
# =============================================================================

def extraer_ca_de_pdb(pdb_file):
    """
    Lee un archivo PDB y devuelve un DataFrame con los Cα de todas las cadenas.

    Columnas: Residuo (nombre), Cadena (ID), Residuo_ID (número), X, Y, Z

    Nota: si el PDB tiene múltiples cadenas (tetrámero, dímero, etc.),
    todos los Cα se incluyen en el mismo DataFrame con su cadena identificada.
    Esto permite calcular contactos tanto intra- como inter-cadena.
    """
    from Bio.PDB import PDBParser
    parser    = PDBParser(QUIET=True)
    structure = parser.get_structure("proteina", pdb_file)

    datos = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):
                    ca = residue["CA"]
                    datos.append({
                        "Residuo"    : residue.resname,
                        "Cadena"     : chain.id,
                        "Residuo_ID" : residue.id[1],
                        "X"          : ca.coord[0],
                        "Y"          : ca.coord[1],
                        "Z"          : ca.coord[2],
                    })
    df = pd.DataFrame(datos)
    print(f"  {pdb_file}: {len(df)} carbonos alfa en {df['Cadena'].nunique()} cadena(s)")
    return df


def calcular_contactos_desde_pdb(pdb_file, output_file,
                                  distancia_umbral=8.0, offset=4):
    """
    Calcula contactos entre todos los pares de Cα en una estructura PDB estática.

    Un par (i, j) se considera contacto si:
        - distancia euclidiana < distancia_umbral (en Å)
        - |posición_global_i - posición_global_j| >= offset
          (para excluir vecinos inmediatos en secuencia)

    Como la estructura es estática, TotalFrac = 1.0 para todos los contactos.
    Esto equivale a asumir que el cristal representa el estado nativo de forma
    perfecta. Si preferís pesos por distancia, podés modificar esa columna.

    Args:
        pdb_file        : archivo PDB de entrada
        output_file     : archivo de salida con formato cpptraj
        distancia_umbral: en Ångströms (8.0 Å es el estándar para modelos Cα)
        offset          : mínima separación en secuencia para considerar contacto

    Output:
        Archivo .dat con columnas: #Res1  #Res2  TotalFrac  Contacts
        Los números de residuo son GLOBALES (no reinician por cadena).
    """
    df = extraer_ca_de_pdb(pdb_file)
    n  = len(df)

    contactos = []
    for i, j in combinations(range(n), 2):
        if j - i < offset:
            continue  # ignorar vecinos cercanos en secuencia

        xi = df.iloc[i][["X","Y","Z"]].values
        xj = df.iloc[j][["X","Y","Z"]].values
        dist = np.linalg.norm(xi - xj)

        if dist < distancia_umbral:
            contactos.append({
                "#Res1"    : i + 1,   # base 1 (convención cpptraj/PDB)
                "#Res2"    : j + 1,
                "TotalFrac": 1.0,     # estructura estática → siempre en contacto
                "Contacts" : 1,
            })

    df_out = pd.DataFrame(contactos)
    df_out.to_csv(output_file, sep=" ", index=False)
    print(f"  → {len(df_out)} contactos guardados en {output_file}")
    return df_out


def ruta_A(pdb_abierto=PDB_ABIERTO, pdb_cerrado=PDB_CERRADO,
           out_abierto=OUT_ABIERTO,  out_cerrado=OUT_CERRADO,
           distancia=DISTANCIA_CONTACTO_A, offset=OFFSET_SECUENCIA):
    """
    Ejecuta la RUTA A completa: PDB → contactos.

    Llama a esta función si partís de archivos PDB cristalográficos.
    """
    print("=== RUTA A: Contactos desde PDB cristalográfico ===")
    print(f"Procesando estado abierto: {pdb_abierto}")
    calcular_contactos_desde_pdb(pdb_abierto, out_abierto,  distancia, offset)
    print(f"Procesando estado cerrado: {pdb_cerrado}")
    calcular_contactos_desde_pdb(pdb_cerrado, out_cerrado, distancia, offset)
    print("RUTA A completada.\n")


# =============================================================================
# RUTA B — Contactos desde trayectoria all-atom (cpptraj)
# =============================================================================

def generar_script_cpptraj_contactos(parm_file, traj_file, output_dat, script_out,
                                      distancia=8, residuos_CA="@CA"):
    """
    Genera el script de cpptraj para calcular contactos nativos sobre
    una trayectoria all-atom.

    El script resultante debe correrse en la terminal con:
        cpptraj -i <script_out>

    Qué hace cpptraj con 'nativecontacts':
        Para cada frame, identifica qué pares de Cα están a menos de
        'distancia' Å y calcula la fracción de frames (TotalFrac) en
        que cada par está en contacto sobre toda la trayectoria.

    Args:
        parm_file   : topología AMBER (.parm7)
        traj_file   : trayectoria (.nc). Puede ser un glob pattern: "sim*.nc"
        output_dat  : nombre del archivo de salida de cpptraj
        script_out  : nombre del script .cpptraj a generar
        distancia   : umbral en Å (por defecto 8, estándar para Cα)
        residuos_CA : máscara de átomos en cpptraj (por defecto solo Cα)
    """
    with open(script_out, "w") as f:
        f.write(f"parm {parm_file}\n")
        f.write(f"trajin {traj_file}\n")
        f.write(f"nativecontacts {residuos_CA} resout {output_dat} distance {distancia}\n")
        f.write("go\nquit\n")
    print(f"  Script cpptraj generado: {script_out}")
    print(f"  Corré: cpptraj -i {script_out}")


def generar_script_cpptraj_extraer_ca(parm_full, traj_glob,
                                       n_residuos_proteina,
                                       parm_ca_out, pdb_frame_out,
                                       frame_numero=777, script_out="preparar_ca.cpptraj"):
    """
    Genera el script de cpptraj para:
        1. Extraer solo los Cα de una trayectoria all-atom (strip)
        2. Guardar el .parm7 y .nc reducidos (solo Cα) → usados después
        3. Guardar un frame específico como PDB → punto de partida de la CG

    Este paso es necesario antes de calcular distancias (RUTA B paso 2)
    porque cpptraj necesita una topología compatible con los índices
    de residuo que aparecen en el .dat de contactos.

    Args:
        parm_full          : topología all-atom completa
        traj_glob          : patrón de trayectoria (ej: "sim*.nc" o "[1-9]00ns.nc")
        n_residuos_proteina: número de residuos a conservar (excluye solvente, etc.)
        parm_ca_out        : nombre del .parm7 reducido a solo Cα
        pdb_frame_out      : nombre del PDB del frame extraído
        frame_numero       : qué frame guardar como referencia (por defecto 777)
        script_out         : nombre del script a generar
    """
    with open(script_out, "w") as f:
        f.write(f"parm {parm_full}\n")
        f.write(f"trajin {traj_glob}\n")
        f.write(f"strip !:(1-{n_residuos_proteina}) parmout {parm_ca_out}\n")
        f.write(f"trajout {pdb_frame_out} start {frame_numero} stop {frame_numero}\n")
        f.write("go\nquit\n")
    print(f"  Script cpptraj generado: {script_out}")
    print(f"  Corré: cpptraj -i {script_out}")


def ruta_B(parm_abierto=PARM_ABIERTO,  traj_abierto=TRAJ_ABIERTO,
           parm_cerrado=PARM_CERRADO,  traj_cerrado=TRAJ_CERRADO,
           out_abierto=OUT_ABIERTO,     out_cerrado=OUT_CERRADO,
           distancia=DISTANCIA_CONTACTO_B):
    """
    Ejecuta la RUTA B: genera los scripts de cpptraj para calcular contactos.

    Después de correr esta función, ejecutá en la terminal:
        cpptraj -i contactos_abierto.cpptraj
        cpptraj -i contactos_cerrado.cpptraj

    Los archivos OUT_ABIERTO y OUT_CERRADO resultantes son los inputs del Paso 2.
    """
    print("=== RUTA B: Contactos desde trayectoria all-atom (cpptraj) ===")
    generar_script_cpptraj_contactos(
        parm_abierto, traj_abierto, out_abierto,
        script_out="contactos_abierto.cpptraj", distancia=distancia
    )
    generar_script_cpptraj_contactos(
        parm_cerrado, traj_cerrado, out_cerrado,
        script_out="contactos_cerrado.cpptraj", distancia=distancia
    )
    print("RUTA B: scripts generados. Corré cpptraj y luego seguí con el Paso 2.\n")


# =============================================================================
# EJECUCIÓN — descomentá la ruta que corresponda
# =============================================================================

if __name__ == "__main__":

    # ── RUTA A: si tenés solo archivos PDB ──────────────────────────────────
    ruta_A()

    # ── RUTA B: si tenés trayectorias all-atom ──────────────────────────────
    # ruta_B()

    # ── Nota: si querés correr RUTA B, primero generá los Cα con cpptraj:
    # generar_script_cpptraj_extraer_ca(
    #     parm_full           = "topologia_full.parm7",
    #     traj_glob           = "[1-9]00ns.nc",
    #     n_residuos_proteina = 1140,          # ajustar a tu sistema
    #     parm_ca_out         = "alfacarbonsAbierto.parm7",
    #     pdb_frame_out       = "single_frame.pdb",
    #     frame_numero        = 777,
    #     script_out          = "preparar_ca.cpptraj"
    # )
