"""
=============================================================================
PASO 2 — FILTRAR CONTACTOS Y OBTENER LAS LISTAS PARA LA SIMULACIÓN CG
=============================================================================

OBJETIVO:
    Tomar los contactos crudos del Paso 1 y producir las dos listas
    que necesita la simulación CG:

    1. contactos_en_comun_openmm.dat  → contactos compartidos por ambos estados
                                        (fuerzas de Gō siempre activas)

    2. contactos_abierto_openmm.dat   → contactos exclusivos del estado abierto
                                        (fuerzas activas o apagadas según variante)

FLUJO DE ESTE SCRIPT:
    [res_contactsAbierto.dat]  ──┐
                                 ├─→ filtrar (TotalFrac > umbral)
    [res_contactsCerrado.dat]  ──┘
            │
            ↓
    [filtrad03abierto.dat]  ──┐
                               ├─→ concatenar → [concat.dat]
    [filtrad03cerrado.dat]  ──┘        │
                                       ↓
                              generar scripts cpptraj
                              para medir distancias en
                              cada estado → correr cpptraj
                                       │
                                       ↓
                              filtrar por ratio de distancias
                              (abierto/cerrado > 1.5)
                                       │
                              ┌────────┴────────────┐
                              ↓                     ↓
                   [contactos_en_comun]    [contactos_exclusivos_abierto]
                              │                     │
                              ↓                     ↓
                   agregar formato OpenMM (columnas de 1s)
                              │                     │
                              ↓                     ↓
              contactos_en_comun_openmm.dat   contactos_abierto_openmm.dat

INPUTS:
    res_contactsAbierto.dat      (del Paso 1)
    res_contactsCerrado.dat      (del Paso 1)
    alfacarbonsAbierto.parm7/.nc (RUTA B) o ca_onlyAbierto.pdb (RUTA A)
    alfacarbonsCerrado.parm7/.nc (RUTA B) o ca_onlyCerrado.pdb (RUTA A)

OUTPUTS:
    contactos_en_comun_openmm.dat   → input directo de la simulación CG
    contactos_abierto_openmm.dat    → input directo de la simulación CG

=============================================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# PARÁMETROS — editá estos valores para tu sistema
# =============================================================================

# Inputs del Paso 1
CONTACTOS_ABIERTO_RAW  = "res_contactsAbierto.dat"
CONTACTOS_CERRADO_RAW  = "res_contactsCerrado.dat"

# Archivos de estructura Cα (para calcular distancias nativas)
# RUTA A: archivos PDB con solo Cα
CA_PDB_ABIERTO  = "1z98_assembly.pdb"
CA_PDB_CERRADO  = "2b5f.pdb"
# RUTA B: topología + trayectoria solo-Cα
CA_PARM_ABIERTO = "alfacarbonsAbierto.parm7"
CA_TRAJ_ABIERTO = "alfacarbonsAbierto.nc"
CA_PARM_CERRADO = "alfacarbonsCerrado.parm7"
CA_TRAJ_CERRADO = "alfacarbonsCerrado.nc"

# Umbral de frecuencia (TotalFrac > UMBRAL para conservar un contacto)
# 0.3 significa: el contacto debe estar formado en al menos 30% de los frames
# Para RUTA A (PDB estático) TotalFrac=1 siempre, así que cualquier umbral < 1 los conserva todos
UMBRAL_TOTALFRAC = 0.3

# Ratio de distancias para identificar contactos "exclusivos del abierto":
# un contacto es "exclusivo del abierto" si dist_abierto > RATIO * dist_cerrado
# Interpretación física: en el estado cerrado ese par está mucho más cerca
# que en el abierto → es un contacto que "se rompe" al abrir
RATIO_DISTANCIAS = 1.5

# Outputs intermedios
OUT_FILTRADO_ABIERTO = "filtrad03_abierto.dat"
OUT_FILTRADO_CERRADO = "filtrad03_cerrado.dat"
OUT_CONCAT           = "contactos_concat.dat"
OUT_DISTANCIAS_SCRIPT_ABIERTO = "distancias_abierto.cpptraj"
OUT_DISTANCIAS_SCRIPT_CERRADO = "distancias_cerrado.cpptraj"
OUT_DIST_ABIERTO_RAW = "contactDistanceabierto.dat"
OUT_DIST_CERRADO_RAW = "contactDistancecerrado.dat"

# Outputs finales (inputs de la simulación)
OUT_COMUN_OPENMM    = "contactos_en_comun_openmm.dat"
OUT_ABIERTO_OPENMM  = "contactos_abierto_openmm.dat"


# =============================================================================
# PASO 2.1 — Filtrar por TotalFrac
# =============================================================================

def filtrar_por_totalfrac(input_file, output_file, umbral=0.3):
    """
    Conserva solo los contactos cuya TotalFrac supera el umbral.

    TotalFrac es la fracción de frames de la trayectoria en que el par
    de residuos estaba en contacto. Un umbral de 0.3 significa que el
    contacto es "persistente" (formado al menos 30% del tiempo).

    Para estructuras estáticas (RUTA A), TotalFrac = 1.0 siempre,
    así que todos los contactos pasan el filtro con cualquier umbral < 1.

    Args:
        input_file  : archivo .dat con columnas #Res1 #Res2 TotalFrac Contacts
        output_file : archivo filtrado
        umbral      : valor mínimo de TotalFrac (por defecto 0.3)
    """
    df = pd.read_csv(input_file, sep=r'\s+')
    df_filtrado = df[df['TotalFrac'] > umbral].copy()
    df_filtrado.to_csv(output_file, sep=' ', index=False)
    print(f"  {input_file}: {len(df)} contactos → {len(df_filtrado)} con TotalFrac > {umbral}")
    return df_filtrado


# =============================================================================
# PASO 2.2 — Intersección: contactos en común entre abierto y cerrado
# =============================================================================

def obtener_contactos_en_comun(file_abierto, file_cerrado, output_file):
    """
    Encuentra los pares (Res1, Res2) que aparecen en AMBOS estados.

    Estos contactos "en común" son los que mantienen la estructura global
    independientemente de la conformación. Se usan como fuerzas de Gō
    siempre activas en la simulación CG (ksb = 1).

    La intersección se hace por merge en (#Res1, #Res2), por lo que
    ambos archivos deben tener la misma numeración de residuos.

    Args:
        file_abierto : contactos filtrados del estado abierto
        file_cerrado : contactos filtrados del estado cerrado
        output_file  : contactos presentes en ambos estados
    """
    df_ab = pd.read_csv(file_abierto, sep=r'\s+')
    df_ce = pd.read_csv(file_cerrado, sep=r'\s+')

    # Renombrar para el merge (los archivos tienen '#Res1', '#Res2')
    comun = pd.merge(df_ab, df_ce, on=['#Res1', '#Res2'])

    # Conservar solo columnas del abierto (TotalFrac_x, Contacts_x)
    comun = comun[['#Res1', '#Res2', 'TotalFrac_x', 'Contacts_x']]
    comun.columns = ['#Res1', '#Res2', 'TotalFrac', 'Contacts']

    comun.to_csv(output_file, sep=' ', index=False)
    print(f"  Contactos en común: {len(comun)} pares → {output_file}")
    return comun


# =============================================================================
# PASO 2.3 — Concatenar ambas listas (para medir distancias de todos)
# =============================================================================

def concatenar_sin_duplicados(file1, file2, output_file):
    """
    Une las listas de contactos de ambos estados eliminando duplicados.

    El resultado es la lista completa de todos los pares que son contacto
    en al menos uno de los dos estados. Se usa para generar el script de
    cpptraj que medirá las distancias de cada par en cada estado.

    Si un par aparece en ambos archivos, se conserva solo una vez
    (la primera aparición, del archivo 1).

    Args:
        file1, file2 : archivos filtrados de cada estado
        output_file  : lista unificada sin duplicados
    """
    df1 = pd.read_csv(file1, sep=' ')
    df2 = pd.read_csv(file2, sep=' ')
    concat = pd.concat([df1, df2], ignore_index=True)
    concat.drop_duplicates(subset=['#Res1', '#Res2'], keep='first', inplace=True)
    concat.to_csv(output_file, index=False, sep=' ', header=False)
    print(f"  Lista concatenada: {len(concat)} pares únicos → {output_file}")
    return concat


# =============================================================================
# PASO 2.4 — Medir distancias en cada estado
# =============================================================================

def generar_script_distancias_cpptraj(input_contactos, parm_file, traj_file,
                                       output_dist, script_out, estado):
    """
    [RUTA B] Genera un script de cpptraj que mide la distancia Cα-Cα
    de cada par de contactos en la trayectoria all-atom.

    Esto permite comparar qué tan separados están los residuos en
    el estado abierto vs. cerrado, para identificar cuáles son
    contactos "exclusivos" de un estado.

    Output de cpptraj: una columna por par, un valor por frame.
    Este script lee el ÚLTIMO frame como valor representativo.

    Corré con: cpptraj -i <script_out>
    """
    with open(script_out, 'w') as f:
        f.write(f"parm {parm_file}\n")
        f.write(f"trajin {traj_file} 777 777\n")   # frame 777 como representativo
        df = pd.read_csv(input_contactos, sep=r'\s+', header=None)
        for _, row in df.iterrows():
            res1, res2 = int(row.iloc[0]), int(row.iloc[1])
            f.write(f"distance :{res1} :{res2} out {output_dist}\n")
    print(f"  Script cpptraj generado: {script_out} → corré: cpptraj -i {script_out}")


def calcular_distancias_desde_pdb(input_contactos, pdb_ca, output_file):
    """
    [RUTA A] Calcula la distancia Cα-Cα de cada par de contactos
    directamente desde el archivo PDB estático.

    Equivalente al script de cpptraj pero sin necesitar trayectoria.
    Lee las coordenadas del PDB y calcula distancias euclidianas.

    Output: una distancia por línea, en el mismo orden que input_contactos.
    """
    import mdtraj as md

    coord = md.load_pdb(pdb_ca)
    df = pd.read_csv(input_contactos, sep=r'\s+', header=None)

    distancias = []
    for _, row in df.iterrows():
        i = int(row.iloc[0]) - 1   # base 0 para MDTraj
        j = int(row.iloc[1]) - 1
        d = md.compute_distances(coord, np.array([[i, j]]))[0][0] * 10  # nm → Å
        distancias.append(d)

    with open(output_file, 'w') as f:
        # Formato que espera transformar_distancias(): 2 líneas (encabezado + datos)
        f.write("frame\n")
        f.write(" ".join(f"{d:.4f}" for d in distancias) + "\n")

    print(f"  {len(distancias)} distancias calculadas → {output_file}")
    return distancias


# =============================================================================
# PASO 2.5 — Transformar y combinar distancias
# =============================================================================

def transformar_distancias(input_file, output_file):
    """
    Convierte el output de cpptraj (o calcular_distancias_desde_pdb)
    de formato horizontal (una fila con todas las distancias) a
    formato vertical (una distancia por línea).

    Esto facilita la comparación uno a uno entre estados.

    Input:  frame1  d1  d2  d3 ...
    Output: d1\n
            d2\n
            d3\n
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()
    # La segunda fila contiene los valores (la primera es encabezado de cpptraj)
    distancias = lines[1].strip().split()
    with open(output_file, 'w') as f:
        for d in distancias:
            f.write(d + '\n')
    print(f"  Distancias transformadas → {output_file} ({len(distancias)} valores)")


def combinar_distancias(file_abierto_transf, file_cerrado_transf, output_file):
    """
    Combina las distancias de ambos estados en un archivo de dos columnas:
        col1 = distancia en estado abierto
        col2 = distancia en estado cerrado

    Una fila por par de contacto, en el mismo orden que la lista concatenada.
    Este archivo se usa para visualizar y filtrar por ratio abierto/cerrado.
    """
    with open(file_abierto_transf) as f1, open(file_cerrado_transf) as f2:
        d_ab = f1.readlines()
        d_ce = f2.readlines()

    if len(d_ab) != len(d_ce):
        raise ValueError(
            f"Número de distancias no coincide: "
            f"abierto={len(d_ab)}, cerrado={len(d_ce)}"
        )

    with open(output_file, 'w') as f:
        for a, c in zip(d_ab, d_ce):
            f.write(f"{a.strip()}\t{c.strip()}\n")
    print(f"  Distancias combinadas → {output_file}")


def visualizar_distancias(distancias_combinadas, ratio=1.5):
    """
    Scatter plot de distancia_abierto vs distancia_cerrado para cada par.

    Los puntos por encima de la línea y = ratio*x corresponden a
    contactos donde la distancia en el abierto es mucho mayor que en
    el cerrado → son candidatos a "contactos exclusivos del cerrado"
    (o dicho de otra forma, son los que SE FORMAN al cerrarse).

    OJO con los ejes: el script original grafica cerrado en X y abierto en Y.
    Un punto sobre y = 1.5x significa dist_abierto > 1.5 * dist_cerrado,
    es decir: en el abierto ese par está más separado → es un contacto
    que existe en el cerrado pero no en el abierto.
    """
    col_ab, col_ce = [], []
    with open(distancias_combinadas) as f:
        for line in f:
            a, c = line.strip().split()
            col_ab.append(float(a))
            col_ce.append(float(c))

    x = np.linspace(min(col_ab + col_ce), max(col_ab + col_ce), 100)

    plt.figure(figsize=(8, 6))
    plt.scatter(col_ce, col_ab, c='steelblue', s=15, alpha=0.6, label='Contactos')
    plt.plot(x, ratio * x, 'g--', label=f'y = {ratio}x (umbral)')
    plt.plot(x, x,          'r-',  label='y = x')
    plt.xlabel('Distancia en estado cerrado (Å)')
    plt.ylabel('Distancia en estado abierto (Å)')
    plt.title('Comparación de distancias por estado conformacional')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    return col_ab, col_ce


# =============================================================================
# PASO 2.6 — Filtrar contactos exclusivos del abierto por ratio de distancias
# =============================================================================

def filtrar_por_ratio_distancias(contactos_concat, dist_abierto_transf,
                                  dist_cerrado_transf, output_exclusivos_abierto,
                                  ratio=1.5):
    """
    Identifica los contactos donde dist_abierto > ratio * dist_cerrado.

    Interpretación:
        Si en el estado abierto dos residuos están mucho más lejos que
        en el cerrado, ese contacto es "exclusivo del cerrado" (o sea,
        se FORMA al cerrar la proteína).

        Estos son los contactos que en la variante "apagada" se ponen
        con ksb=0: si los apagás, el sistema no tiene incentivo energético
        para cerrar y permanece en el estado abierto.

    Nota: el script original usa la condición invertida (dist_cerrado > ratio * dist_abierto)
    para definir los "exclusivos del abierto". Verificá cuál es la convención
    de tu sistema mirando el scatter plot del paso anterior.

    Args:
        contactos_concat          : lista completa de pares (sin encabezado)
        dist_abierto_transf       : distancias verticales del estado abierto
        dist_cerrado_transf       : distancias verticales del estado cerrado
        output_exclusivos_abierto : archivo de salida
        ratio                     : umbral (por defecto 1.5)
    """
    df = pd.read_csv(contactos_concat, sep=r'\s+', header=None)

    with open(dist_abierto_transf) as f:
        d_ab = [float(l.strip()) for l in f]
    with open(dist_cerrado_transf) as f:
        d_ce = [float(l.strip()) for l in f]

    # Condición: dist_abierto > ratio * dist_cerrado
    # → el par está más lejos en el abierto → es exclusivo del cerrado
    # → si lo "apagás" el sistema pierde la fuerza que lo cierra
    mask = [a > ratio * c for a, c in zip(d_ab, d_ce)]
    df_excl = df[mask].copy()

    df_excl.to_csv(output_exclusivos_abierto, sep=' ', index=False, header=False)
    print(f"  Contactos filtrados (ratio > {ratio}): {df_excl.sum()} pares → {output_exclusivos_abierto}")
    return df_excl


# =============================================================================
# PASO 2.7 — Convertir al formato que espera la simulación OpenMM
# =============================================================================

def convertir_a_formato_openmm(input_file, output_file):
    """
    Agrega dos columnas de "1" en las posiciones 0 y 2, que es el
    formato que espera generar_lista_contactos() en la simulación.

    Formato entrada:  Res1  Res2  TotalFrac  Contacts
    Formato salida:   1     Res1  1          Res2

    Este formato particular viene del script original donde la función
    generar_lista_contactos() lee las columnas 1 y 3 (base 0).

    Si en el futuro modificás generar_lista_contactos(), podés
    simplificar este paso.
    """
    df = pd.read_csv(input_file, sep=r'\s+', header=None)

    # Insertar columnas de 1s en posición 0 y 2
    df.insert(0, 'c0', 1)
    df.insert(2, 'c2', 1)

    # Eliminar columnas sobrantes (TotalFrac y Contacts originales)
    df = df.iloc[:, :4]

    df.to_csv(output_file, sep=' ', index=False, header=False)
    print(f"  Formato OpenMM → {output_file} ({len(df)} contactos)")


# =============================================================================
# PIPELINE COMPLETO
# =============================================================================

def pipeline_completo(modo="A"):
    """
    Ejecuta el pipeline completo del Paso 2.

    Args:
        modo : "A" para RUTA A (desde PDB), "B" para RUTA B (desde trayectoria)

    Al final produce:
        contactos_en_comun_openmm.dat
        contactos_abierto_openmm.dat
    """
    print(f"\n{'='*60}")
    print(f"PASO 2 — Filtrado y preparación de contactos (MODO {modo})")
    print(f"{'='*60}\n")

    # 2.1 — Filtrar por TotalFrac
    print("── 2.1 Filtrando por TotalFrac ──")
    filtrar_por_totalfrac(CONTACTOS_ABIERTO_RAW, OUT_FILTRADO_ABIERTO, UMBRAL_TOTALFRAC)
    filtrar_por_totalfrac(CONTACTOS_CERRADO_RAW, OUT_FILTRADO_CERRADO, UMBRAL_TOTALFRAC)

    # 2.2 — Contactos en común
    print("\n── 2.2 Intersección (contactos en común) ──")
    obtener_contactos_en_comun(OUT_FILTRADO_ABIERTO, OUT_FILTRADO_CERRADO,
                                "contactos_en_comun.dat")

    # 2.3 — Concatenar
    print("\n── 2.3 Concatenando listas ──")
    concatenar_sin_duplicados(OUT_FILTRADO_ABIERTO, OUT_FILTRADO_CERRADO, OUT_CONCAT)

    # 2.4 — Medir distancias en cada estado
    print(f"\n── 2.4 Midiendo distancias (MODO {modo}) ──")
    if modo == "A":
        # RUTA A: calcular desde PDB estático
        calcular_distancias_desde_pdb(OUT_CONCAT, CA_PDB_ABIERTO, OUT_DIST_ABIERTO_RAW)
        calcular_distancias_desde_pdb(OUT_CONCAT, CA_PDB_CERRADO, OUT_DIST_CERRADO_RAW)
    else:
        # RUTA B: generar scripts de cpptraj (hay que correrlos manualmente)
        generar_script_distancias_cpptraj(
            OUT_CONCAT, CA_PARM_ABIERTO, CA_TRAJ_ABIERTO,
            OUT_DIST_ABIERTO_RAW, OUT_DISTANCIAS_SCRIPT_ABIERTO, "abierto"
        )
        generar_script_distancias_cpptraj(
            OUT_CONCAT, CA_PARM_CERRADO, CA_TRAJ_CERRADO,
            OUT_DIST_CERRADO_RAW, OUT_DISTANCIAS_SCRIPT_CERRADO, "cerrado"
        )
        print("  *** DETENÉ aquí y corré los scripts de cpptraj ***")
        print("  Luego volvé a correr a partir de transformar_distancias()")
        return   # salir hasta que el usuario corra cpptraj

    # 2.5 — Transformar y combinar
    print("\n── 2.5 Transformando y combinando distancias ──")
    DIST_AB_TRANSF = "dist_abierto_transf.dat"
    DIST_CE_TRANSF = "dist_cerrado_transf.dat"
    DIST_COMBINADA = "distancias_combinadas.dat"

    transformar_distancias(OUT_DIST_ABIERTO_RAW, DIST_AB_TRANSF)
    transformar_distancias(OUT_DIST_CERRADO_RAW, DIST_CE_TRANSF)
    combinar_distancias(DIST_AB_TRANSF, DIST_CE_TRANSF, DIST_COMBINADA)

    # Visualización
    col_ab, col_ce = visualizar_distancias(DIST_COMBINADA, ratio=RATIO_DISTANCIAS)

    # 2.6 — Filtrar por ratio
    print("\n── 2.6 Filtrando por ratio de distancias ──")
    CONTACTOS_EXCL_ABIERTO = "contactos_exclusivos_abierto.dat"
    filtrar_por_ratio_distancias(
        OUT_CONCAT, DIST_AB_TRANSF, DIST_CE_TRANSF,
        CONTACTOS_EXCL_ABIERTO, ratio=RATIO_DISTANCIAS
    )

    # 2.7 — Formato OpenMM
    print("\n── 2.7 Convirtiendo a formato OpenMM ──")
    convertir_a_formato_openmm("contactos_en_comun.dat", OUT_COMUN_OPENMM)
    convertir_a_formato_openmm(CONTACTOS_EXCL_ABIERTO,  OUT_ABIERTO_OPENMM)

    print(f"\n{'='*60}")
    print("PASO 2 COMPLETADO. Archivos listos para la simulación:")
    print(f"  → {OUT_COMUN_OPENMM}")
    print(f"  → {OUT_ABIERTO_OPENMM}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    # Elegí "A" si partís de PDB, "B" si tenés trayectorias all-atom
    pipeline_completo(modo="A")
