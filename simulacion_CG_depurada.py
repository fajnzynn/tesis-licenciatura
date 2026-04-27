"""
=============================================================================
SIMULACIÓN DE DINÁMICA MOLECULAR A GRANO GRUESO (CG) CON OPENMM
Modelo tipo Gō basado en contactos — un carbono alfa (Cα) por residuo
=============================================================================

DESCRIPCIÓN GENERAL DEL MODELO:
    - Cada residuo está representado por una única partícula ubicada en su Cα.
    - Las fuerzas que mantienen la estructura son:
        1. Enlaces Cα-Cα consecutivos (HarmonicBondForce)
        2. Ángulos entre tres Cα consecutivos (HarmonicAngleForce)
        3. Diedros entre cuatro Cα consecutivos (CustomTorsionForce)
        4. Contactos estructura-based (fuerzas de Gō, CustomBondForce)
        5. Exclusión estérica entre pares no-nativos (repulsión suave)

INPUTS REQUERIDOS:
    - single_frame.pdb       : estructura all-atom de entrada
    - ca_onlyCerrado.pdb     : estructura de referencia cerrada (solo Cα)
    - ca_onlyAbierto.pdb     : estructura de referencia abierta (solo Cα)
    - contactos_en_comun03_openmm.dat : contactos compartidos abierto+cerrado
    - contactos1_5_abierto.dat        : contactos exclusivos del estado abierto
    - cg.xml                 : forcefield coarse-grained (define partículas CG)

OUTPUTS:
    - out-fric2-N.dcd   : trayectorias de cada réplica
    - lastconf.pdb      : última conformación guardada
    - data.an           : energía potencial y temperatura en función del tiempo
"""

# =============================================================================
# SECCIÓN 1: IMPORTS
# =============================================================================

from openmm.app import *
from openmm import *
from openmm.app import PDBFile
from pdbfixer import PDBFixer
import numpy as np
import mdtraj as md
from openmm import CustomBondForce
from openmm.unit import *
from openmm import unit
from sys import stdout
import random as rd
import matplotlib.pyplot as plt
import pandas as pd

# Nota: si usás una versión antigua de OpenMM (< 7.6), los imports son:
# from simtk.openmm.app import * / from simtk.openmm import *

# =============================================================================
# SECCIÓN 2: ARCHIVOS DE ENTRADA Y SALIDA
# =============================================================================

# Forcefield coarse-grained: define la masa y tipo de partícula para cada residuo CG
FORCEFIELD_FILE = "cg.xml"

# Estructura all-atom de entrada (puede tener heteroátomos, agua, etc.)
INPUT_PDB       = "1z98_assembly.pdb"

# Archivos intermedios generados en el preprocesamiento
PROTEIN_PDB     = "protein_abierto.pdb"       # Estructura limpia (sin heteroátomos)
CA_ONLY_PDB     = "ca_only_abierto.pdb"  # Solo carbonos alfa — input de la simulación

# Estructuras de referencia para calcular distancias nativas (r_ijN)
ESTRUCTURA_CERRADA = "ca_only_cerrado.pdb"
ESTRUCTURA_ABIERTA = "ca_only_abierto.pdb"

# Archivos de contactos (generados por el pipeline de análisis)
CONTACTOS_EN_COMUN  = "contactos_en_comun_openmm.dat"  # Compartidos abierto+cerrado
CONTACTOS_ABIERTO   = "contactos_abierto_openmm.dat"          # Exclusivos del abierto

forcefield = ForceField(FORCEFIELD_FILE)

# =============================================================================
# SECCIÓN 3: PREPROCESAMIENTO — extraer carbonos alfa desde la estructura all-atom
# =============================================================================

def preparar_estructura(input_pdb, output_cleaned, output_ca):
    """
    1. Limpia la estructura all-atom: elimina agua y heteroátomos.
    2. Extrae solo los carbonos alfa y los guarda en un PDB separado.

    Args:
        input_pdb      : archivo PDB all-atom original
        output_cleaned : PDB limpio (sin heteroátomos)
        output_ca      : PDB con solo los Cα
    """
    # Limpiar estructura con PDBFixer
    fixer = PDBFixer(input_pdb)
    fixer.removeHeterogens(keepWater=False)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_cleaned, 'w'))

    # Extraer solo los átomos Cα con MDTraj
    pdb = md.load(output_cleaned)
    indices_ca = [atom.index for atom in pdb.topology.atoms if atom.name == "CA"]
    solo_ca = pdb.atom_slice(indices_ca)
    solo_ca.save(output_ca)
    print(f"Preprocesamiento completo: {len(indices_ca)} carbonos alfa extraídos.")

preparar_estructura(INPUT_PDB, PROTEIN_PDB, CA_ONLY_PDB)

# =============================================================================
# SECCIÓN 4: LECTURA DE CONTACTOS
# =============================================================================

def generar_lista_contactos(contactFile):
    """
    Lee un archivo de contactos con formato:
        TotalFrac  Res1  Contacts  Res2
    y devuelve una lista de pares (i, j) en índice base 0.

    El archivo tiene una fila de encabezado que se saltea.
    Los números de residuo en el archivo son base 1, se convierten a base 0.
    """
    df = pd.read_csv(contactFile, sep=" ", skiprows=1, header=None)
    df.columns = ['totalfrac', 'i', 'contacts', 'j']
    df = df[['i', 'j']].astype(int)
    # Convertir de base 1 (numeración PDB) a base 0 (índice Python)
    return list(zip(df['i'] - 1, df['j'] - 1))


def generar_lista_parametros(contact_list, pdb_referencia):
    """
    Para cada par de contacto (i, j), calcula la distancia nativa r_ijN
    en la estructura de referencia dada.

    Estas distancias son los parámetros de equilibrio del potencial de Gō.

    Args:
        contact_list   : lista de pares (i, j) en base 0
        pdb_referencia : archivo PDB con solo Cα de la conformación de referencia

    Returns:
        Lista de distancias nativas en nanómetros
    """
    coord = md.load_pdb(pdb_referencia)
    parameter_list = []
    for i, j in contact_list:
        dist = md.compute_distances(coord, np.array([[i, j]]))[0][0]
        parameter_list.append(dist)
    return parameter_list


def visualizar_mapa_contactos(contact_list, parameter_list, titulo="Mapa de contactos"):
    """
    Genera un scatter plot del mapa de contactos coloreado por distancia nativa.
    Útil para verificar visualmente que los contactos tienen sentido estructural.
    """
    x = [par[0] for par in contact_list]
    y = [par[1] for par in contact_list]
    plt.figure()
    sc = plt.scatter(x, y, c=parameter_list, cmap='viridis')
    plt.colorbar(sc, label='Distancia nativa (nm)')
    plt.xlabel('Residuo i')
    plt.ylabel('Residuo j')
    plt.title(titulo)
    plt.show()


# Cargar contactos del estado abierto (exclusivos): se usan con ksb=0 en variante apagada
lista_contactos_abierto   = generar_lista_contactos(CONTACTOS_ABIERTO)
lista_parametros_abierto  = generar_lista_parametros(lista_contactos_abierto, ESTRUCTURA_ABIERTA)

# Cargar contactos en común (abierto + cerrado): siempre activos con ksb=1
lista_contactos_en_comun  = generar_lista_contactos(CONTACTOS_EN_COMUN)
lista_parametros_en_comun = generar_lista_parametros(lista_contactos_en_comun, ESTRUCTURA_CERRADA)

visualizar_mapa_contactos(lista_contactos_abierto,  lista_parametros_abierto,  "Contactos abierto")
visualizar_mapa_contactos(lista_contactos_en_comun, lista_parametros_en_comun, "Contactos en común")

# =============================================================================
# SECCIÓN 5: CÁLCULO DE GEOMETRÍA NATIVA (enlaces, ángulos, diedros)
# =============================================================================

def calcular_geometria_nativa(pdb_ca):
    """
    A partir de la estructura Cα de referencia, extrae:
      - Distancias de enlace entre Cα consecutivos dentro de cada cadena
      - Ángulos entre tres Cα consecutivos dentro de cada cadena
      - Diedros entre cuatro Cα consecutivos dentro de cada cadena

    Los términos de enlace, ángulo y diedro son ARMÓNICOS alrededor de estos
    valores nativos, lo que mantiene la topología de cadena correcta.

    También devuelve chainRanges: lista de tuplas (inicio, fin) con los índices
    globales de cada cadena, necesaria para no conectar el último residuo de
    una cadena con el primero de la siguiente.

    Returns:
        enlaces, angulos, diedros : arrays de MDTraj con valores nativos
        chainRanges               : lista de tuplas (first, last) por cadena
    """
    pdb = md.load_pdb(pdb_ca)
    topology = pdb.topology

    bonds_indices    = []
    angles_indices   = []
    dihedral_indices = []
    chainRanges      = []

    offset = 0  # índice global acumulado
    for chain in topology.chains:
        residuos = list(chain.residues)
        n = len(residuos)
        first = offset
        last  = offset + n
        chainRanges.append((first, last))

        for i in range(first, last - 1):
            bonds_indices.append((i, i+1))
        for i in range(first, last - 2):
            angles_indices.append((i, i+1, i+2))
        for i in range(first, last - 3):
            dihedral_indices.append((i, i+1, i+2, i+3))

        offset += n

    enlaces  = md.compute_distances(pdb, bonds_indices,    periodic=False)
    angulos  = md.compute_angles(   pdb, angles_indices,   periodic=False)
    diedros  = md.compute_dihedrals(pdb, dihedral_indices, periodic=False)

    print(f"Cadenas: {len(chainRanges)}")
    for idx, (f, l) in enumerate(chainRanges):
        print(f"  Cadena {idx}: residuos {f} a {l-1} ({l-f} residuos)")

    return enlaces, angulos, diedros, chainRanges


enlaces, angulos_nat, diedros_nat, chainRanges = calcular_geometria_nativa(CA_ONLY_PDB)

# =============================================================================
# SECCIÓN 6: DEFINICIÓN DE LAS FUERZAS DEL MODELO CG
# =============================================================================

def ca_bonds(system, enlaces, k_bond, chainRanges):
    """
    FUERZA 1: Enlace armónico entre Cα consecutivos.

    Potencial: V = k/2 * (r - r0)^2
    donde r0 es la distancia nativa entre Cα_i y Cα_{i+1}.

    k_bond típico: 20000 kJ/mol/nm² — valor alto para mantener la cadena rígida.
    Solo conecta residuos dentro de la misma cadena (usa chainRanges).
    """
    bond_force = HarmonicBondForce()
    bond_idx = 0
    for first, last in chainRanges:
        for i in range(first, last - 1):
            r0 = float(enlaces[0][bond_idx])  # distancia nativa en nm
            bond_force.addBond(i, i+1, r0, k_bond)
            bond_idx += 1
    return bond_force


def angle_term(system, angulos, k_ang, chainRanges):
    """
    FUERZA 2: Ángulo armónico entre tres Cα consecutivos.

    Potencial: V = k/2 * (theta - theta0)^2
    donde theta0 es el ángulo nativo.

    k_ang típico: 40 kJ/mol/rad² — controla la rigidez local de la cadena.
    Solo aplica dentro de la misma cadena.
    """
    angle_force = HarmonicAngleForce()
    ang_idx = 0
    for first, last in chainRanges:
        for i in range(first, last - 2):
            theta0 = float(angulos[0][ang_idx])  # ángulo nativo en radianes
            angle_force.addAngle(i, i+1, i+2, theta0, k_ang)
            ang_idx += 1
    return angle_force


def dihedral_term(system, diedros, k_dih, chainRanges):
    """
    FUERZA 3: Término diedro con dos armónicos (fundamental + tercer armónico).

    Potencial: V = k*(1 - cos(phi - phi0)) + 0.5*k*(1 - cos(3*(phi - phi0)))
    donde phi0 es el diedro nativo.

    Este potencial tiene un mínimo en phi0 y barreras más pequeñas en ±120°,
    lo que permite transiciones conformacionales controladas.

    k_dih fijo en 1 kJ/mol (codificado en la función, no como parámetro).
    Solo aplica dentro de la misma cadena.
    """
    k_dih_val = 1 * (kilojoule_per_mole / radian)
    dih_force = CustomTorsionForce(
        "k_dih*(1-cos(theta-theta0)) + 0.5*k_dih*(1-cos(3*(theta-theta0)))"
    )
    dih_force.addPerTorsionParameter("k_dih")
    dih_force.addPerTorsionParameter("theta0")

    tor_idx = 0
    for first, last in chainRanges:
        for i in range(first, last - 3):
            phi0 = float(diedros[0][tor_idx])  # diedro nativo en radianes
            dih_force.addTorsion(i, i+1, i+2, i+3, [k_dih_val, phi0])
            tor_idx += 1
    return dih_force


def structure_based_term(contact_list, parameter_list, ksb):
    """
    FUERZA 4: Potencial de Gō (estructura-based) para contactos nativos.

    Potencial tipo 10-12 de Lennard-Jones modificado:
        V = ksb * [5*(r_ijN/r)^12 - 6*(r_ijN/r)^10]

    - Mínimo en r = r_ijN (distancia nativa): V = -ksb
    - Repulsivo a r < r_ijN, atractivo a r > r_ijN con decaimiento suave
    - ksb = 0 desactiva completamente el contacto (usado en variante apagada)
    - ksb = 1 kJ/mol es el valor estándar

    Args:
        contact_list   : lista de pares (i, j)
        parameter_list : lista de distancias nativas r_ijN (nm)
        ksb            : energía de contacto en kJ/mol (0 = apagado, 1 = activo)
    """
    go_force = CustomBondForce(
        f"{ksb}*((5*((r_ijN/r)^12)) - (6*((r_ijN/r)^10)))"
    )
    go_force.addPerBondParameter("ksb")
    go_force.addPerBondParameter("r_ijN")

    for (i, j), r_nat in zip(contact_list, parameter_list):
        go_force.addBond(i, j, [ksb * kilojoule_per_mole, r_nat])
    return go_force


def exclusion_term(contact_list, kex, system):
    """
    FUERZA 5: Repulsión estérica entre pares no-nativos.

    Potencial: V = kex * (0.4/r)^12

    Se aplica a todos los pares (i, j) con j >= i+4 que NO estén en
    la lista de contactos nativos. Evita que residuos no contactados
    se sobrepongan en el espacio, dando volumen excluido a la cadena.

    kex típico: 1 kJ/mol
    El umbral j >= i+4 excluye los vecinos en secuencia (ya cubiertos
    por enlace, ángulo y diedro).
    """
    excl_force = CustomBondForce(f"{kex}*((0.4/r)^12)")
    excl_force.addPerBondParameter("kex")
    contact_set = set(contact_list)
    n = system.getNumParticles()
    for i in range(n):
        for j in range(i + 4, n):
            if (i, j) not in contact_set:
                excl_force.addBond(i, j, [kex * kilojoule_per_mole])
    return excl_force

# =============================================================================
# SECCIÓN 7: VARIABLES COLECTIVAS — parámetro de orden Q y Umbrella Sampling
# =============================================================================

def Qcalc(eps, contact_list, parameter_list, gamma):
    """
    Calcula el parámetro de orden Q: fracción de contactos nativos formados.

    Q = (1/N) * sum_i [ 0.5 * (1 - tanh(gamma*(r_i - 1.2*r_ijN))) ]

    Cada término vale ~1 cuando el contacto está formado (r < 1.2*r_ijN)
    y ~0 cuando está roto. El factor 1.2 da tolerancia al contacto.

    Q es una variable colectiva continua y diferenciable, útil como
    coordenada de reacción para estudiar el plegamiento.

    Args:
        gamma : suavidad de la función (típico: 50 nm^-1)
        eps   : escala de energía (típico: 1)
    """
    N = len(contact_list)
    Q_force = CustomBondForce(
        f"1/{N}*eps*0.5*(1-tanh(gamma*(r-1.2*r_ijN)))"
    )
    Q_force.addGlobalParameter("gamma", gamma)
    Q_force.addPerBondParameter("eps")
    Q_force.addPerBondParameter("r_ijN")
    for (i, j), r_nat in zip(contact_list, parameter_list):
        Q_force.addBond(i, j, [eps, r_nat])
    Q_force.setForceGroup(12)  # grupo separado para leer Q sin recomputar todo
    return Q_force


def umbrella_sampling(kq, q0, eps, contact_list, parameter_list, gamma):
    """
    Aplica un potencial armónico sobre Q para forzar la simulación
    a explorar un valor específico de Q (umbrella sampling).

    V_bias = 0.5 * kq * (Q - q0)^2

    Combinando múltiples ventanas (q0 = 0.0, 0.05, ..., 1.0) y
    postprocesando con WHAM, se obtiene el perfil de energía libre F(Q).

    Args:
        kq : constante de fuerza del sesgo (típico: 4000 kJ/mol)
        q0 : valor objetivo de Q para esta ventana
    """
    Q_cv   = Qcalc(eps, contact_list, parameter_list, gamma)
    bias   = CustomCVForce(f"0.5*{kq}*(q-{q0})^2")
    bias.addGlobalParameter("kq", kq * kilojoule_per_mole)
    bias.addGlobalParameter("q0", q0)
    bias.addCollectiveVariable("q", Q_cv)
    bias.setForceGroup(13)
    return bias

# =============================================================================
# SECCIÓN 8: CONSTRUCCIÓN DEL SISTEMA
# =============================================================================

def construir_sistema(ca_pdb, forcefield, enlaces, angulos_nat, diedros_nat,
                      chainRanges,
                      lista_contactos_abierto,   lista_parametros_abierto,
                      lista_contactos_en_comun,  lista_parametros_en_comun,
                      ksb_abierto=0, ksb_comun=1, k_bond=20000, k_ang=40):
    """
    Ensambla el sistema OpenMM con todas las fuerzas del modelo CG.

    Parámetros importantes:
        ksb_abierto : energía de los contactos exclusivos del abierto.
                      0 = apagados (conformación cerrada favorecida)
                      1 = activos  (conformación abierta favorecida)
        ksb_comun   : energía de los contactos en común. Siempre 1.
        k_bond      : rigidez de enlace (kJ/mol/nm²)
        k_ang       : rigidez de ángulo (kJ/mol/rad²)

    Grupos de fuerza (para leer energías por separado):
        1 = enlaces
        2 = ángulos
        3 = diedros
        4 = contactos Gō
        5 = exclusión estérica
        12 = Q (parámetro de orden)
        13 = sesgo de umbrella
    """
    pdb    = PDBFile(ca_pdb)
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=3 * nanometers
    )
    # Eliminar la fuerza nonbonded por defecto del forcefield CG
    system.removeForce(0)

    # --- Fuerzas de cadena ---
    f_bond = ca_bonds(system, enlaces, k_bond, chainRanges)
    system.addForce(f_bond)
    f_bond.setForceGroup(1)

    f_ang = angle_term(system, angulos_nat, k_ang, chainRanges)
    system.addForce(f_ang)
    f_ang.setForceGroup(2)

    f_dih = dihedral_term(system, diedros_nat, 1, chainRanges)
    system.addForce(f_dih)
    f_dih.setForceGroup(3)

    # --- Contactos nativos (fuerzas de Gō) ---
    # Contactos exclusivos del abierto (ksb_abierto=0 en la variante "apagada")
    f_go_abierto = structure_based_term(
        lista_contactos_abierto, lista_parametros_abierto, ksb_abierto
    )
    system.addForce(f_go_abierto)

    # Contactos compartidos (siempre activos)
    f_go_comun = structure_based_term(
        lista_contactos_en_comun, lista_parametros_en_comun, ksb_comun
    )
    system.addForce(f_go_comun)
    f_go_comun.setForceGroup(4)

    # --- Exclusión estérica entre pares no-nativos ---
    todos_contactos = lista_contactos_abierto + lista_contactos_en_comun
    f_excl = exclusion_term(todos_contactos, 1, system)
    system.addForce(f_excl)
    f_excl.setForceGroup(5)

    print(f"Sistema construido: {system.getNumParticles()} partículas, "
          f"{system.getNumForces()} fuerzas")
    return system, pdb


system, pdb_openmm = construir_sistema(
    CA_ONLY_PDB, forcefield, enlaces, angulos_nat, diedros_nat, chainRanges,
    lista_contactos_abierto,  lista_parametros_abierto,
    lista_contactos_en_comun, lista_parametros_en_comun,
    ksb_abierto=0,   # <-- cambiar a 1 para activar contactos del abierto
    ksb_comun=1
)

# =============================================================================
# SECCIÓN 9: SIMULACIÓN — réplicas a temperatura constante (NVT)
# =============================================================================

def correr_replicas_NVT(system, pdb_openmm, n_replicas=20,
                         temperatura=70, n_pasos=100_000,
                         prefijo_salida="out-fric2"):
    """
    Corre n_replicas independientes a temperatura constante con integrador
    de Langevin.

    Cada réplica:
        - Parte de la última conformación de la réplica anterior (encadenadas)
        - Usa una semilla aleatoria diferente para las velocidades iniciales
        - Guarda trayectoria DCD cada 100 pasos
        - Acumula energía y temperatura en data.an

    Parámetros del integrador Langevin:
        - temperatura   : en Kelvin (unidades reducidas CG, no Kelvin real)
        - fricción      : 1/picosecond (amortiguamiento moderado)
        - paso de tiempo: 0.0005 ps

    GPU (CUDA) con precisión mixta para mayor velocidad.
    """
    import os
    data_file = "data.an"
    open(data_file, 'w').close()  # limpiar archivo de salida

    pdb_actual = pdb_openmm  # punto de partida

    for i in range(n_replicas):
        nombre_dcd = f"{prefijo_salida}-{i}.dcd"
        nombre_txt = "data.txt"

        integrador = LangevinIntegrator(
            temperatura * kelvin,
            1 / picoseconds,       # coeficiente de fricción
            0.0005 * picoseconds   # paso de tiempo
        )

        # Usar GPU si está disponible
        try:
            plataforma  = Platform.getPlatformByName('CUDA')
            propiedades = {'CudaPrecision': 'mixed'}
            simulacion  = Simulation(pdb_actual.topology, system,
                                     integrador, plataforma, propiedades)
        except Exception:
            print("CUDA no disponible, usando CPU")
            simulacion = Simulation(pdb_actual.topology, system, integrador)

        simulacion.context.setPositions(pdb_actual.positions)
        simulacion.minimizeEnergy()

        # Velocidades aleatorias con semilla diferente en cada réplica
        semilla = rd.randint(100_000, 1_000_000)
        simulacion.context.setVelocitiesToTemperature(temperatura, semilla)

        # Reporters: trayectoria + datos termodinámicos
        simulacion.reporters.append(DCDReporter(nombre_dcd, 100, append=False))
        simulacion.reporters.append(PDBReporter("lastconf.pdb", n_pasos))
        simulacion.reporters.append(StateDataReporter(
            nombre_txt, 100, step=True,
            potentialEnergy=True, temperature=True
        ))

        simulacion.step(n_pasos)

        # Acumular datos en data.an (salteando encabezado)
        os.system(f"cat {nombre_txt} | tr ',' ' ' | awk 'NR>1' >> {data_file}")

        # Cargar última conformación como punto de partida de la siguiente réplica
        pdb_actual = PDBFile("lastconf.pdb")
        print(f"Réplica {i+1}/{n_replicas} completada → {nombre_dcd}")


correr_replicas_NVT(system, pdb_openmm, n_replicas=20, temperatura=70, n_pasos=100_000)

# =============================================================================
# SECCIÓN 10: SIMULATED ANNEALING — enfriamiento gradual desde T_inicial
# =============================================================================

def simulated_annealing(system, pdb_openmm, T_inicial=120, n_ciclos=80,
                         n_pasos=100_000, delta_T=1):
    """
    Realiza un protocolo de recocido simulado: baja la temperatura de forma
    gradual desde T_inicial hasta ~T_final = T_inicial - n_ciclos * delta_T.

    Útil para:
        - Encontrar conformaciones de mínima energía
        - Escapar de mínimos locales antes de correr producción
        - Explorar el espacio conformacional de forma guiada

    Cada ciclo corre n_pasos con temperatura Tf = T_inicial - delta_T * ciclo.
    """
    import os
    data_file = "data.an"
    open(data_file, 'w').close()

    pdb_actual = pdb_openmm

    for i in range(n_ciclos):
        Tf = T_inicial - delta_T * i   # temperatura decrece linealmente
        nombre_dcd = f"annealing-T{Tf:.0f}.dcd"

        integrador = LangevinIntegrator(
            Tf * kelvin, 1 / picoseconds, 0.0005 * picoseconds
        )

        try:
            plataforma  = Platform.getPlatformByName('CUDA')
            propiedades = {'CudaPrecision': 'mixed'}
            simulacion  = Simulation(pdb_actual.topology, system,
                                     integrador, plataforma, propiedades)
        except Exception:
            simulacion = Simulation(pdb_actual.topology, system, integrador)

        simulacion.context.setPositions(pdb_actual.positions)
        semilla = rd.randint(100_000, 1_000_000)
        simulacion.context.setVelocitiesToTemperature(Tf, semilla)
        simulacion.minimizeEnergy()

        simulacion.reporters.append(DCDReporter(nombre_dcd, 1000, append=False))
        simulacion.reporters.append(PDBReporter("lastconf.pdb", n_pasos))
        simulacion.reporters.append(StateDataReporter(
            "data.txt", 1000, step=True, potentialEnergy=True, temperature=True
        ))
        simulacion.step(n_pasos)

        os.system("cat data.txt | tr ',' ' ' | awk 'NR>1' >> data.an")
        pdb_actual = PDBFile("lastconf.pdb")
        print(f"Annealing ciclo {i+1}/{n_ciclos} — T={Tf:.0f} K → {nombre_dcd}")


# simulated_annealing(system, pdb_openmm, T_inicial=120, n_ciclos=80)
# (descomentar para correr)

# =============================================================================
# SECCIÓN 11: UMBRELLA SAMPLING — muestreo a lo largo de Q
# =============================================================================

def correr_umbrella_sampling(system, pdb_openmm, forcefield,
                              contact_list, parameter_list,
                              Q_values, temperaturas, Kq=4000,
                              n_pasos=250_000):
    """
    Recorre ventanas de umbrella sampling variando Q_objetivo (q0).

    Para cada combinación (q0, temperatura):
        1. Construye el sistema con sesgo armónico en Q
        2. Corre n_pasos de Langevin
        3. Guarda trayectoria como "{q0:.4f}_{temperatura:.4f}.dcd"

    Las trayectorias resultantes se analizan con WHAM para obtener F(Q).

    Args:
        Q_values    : lista/array de valores de Q objetivo (ej: np.arange(0.7, 1.0, 0.025))
        temperaturas: lista de temperaturas a explorar (ej: [127, 130, 133])
        Kq          : constante de fuerza del sesgo en kJ/mol
    """
    pdb_ref = PDBFile(CA_ONLY_PDB)

    for q0 in Q_values:
        # Construir sistema fresco con sesgo de umbrella
        sys_umb, pdb_umb = construir_sistema(
            CA_ONLY_PDB, forcefield, enlaces, angulos_nat, diedros_nat, chainRanges,
            lista_contactos_abierto,  lista_parametros_abierto,
            lista_contactos_en_comun, lista_parametros_en_comun,
        )

        # Agregar Q monitor y sesgo
        Q_monitor = Qcalc(1, contact_list, parameter_list, gamma=50)
        sys_umb.addForce(Q_monitor)

        sesgo = umbrella_sampling(Kq, q0, 1, contact_list, parameter_list, gamma=50)
        sys_umb.addForce(sesgo)

        for T in temperaturas:
            nombre_dcd = f"{q0:.4f}_{T:.4f}.dcd"
            integrador = LangevinIntegrator(
                T * kelvin, 0.5 / picosecond, 0.0005 * picoseconds
            )
            simulacion = Simulation(pdb_umb.topology, sys_umb, integrador)
            simulacion.context.setPositions(pdb_umb.positions)
            simulacion.minimizeEnergy()
            simulacion.reporters.append(DCDReporter(nombre_dcd, 1000, append=False))
            simulacion.step(n_pasos)
            print(f"Umbrella q0={q0:.3f} T={T} K → {nombre_dcd}")


# Ejemplo de uso:
# correr_umbrella_sampling(
#     system, pdb_openmm, forcefield,
#     contact_list  = lista_contactos_en_comun,
#     parameter_list= lista_parametros_en_comun,
#     Q_values      = np.arange(0.7, 1.0, 0.025),
#     temperaturas  = [127, 130, 133],
#     Kq=4000, n_pasos=250_000
# )

# =============================================================================
# SECCIÓN 12: ANÁLISIS POST-SIMULACIÓN — energías y Q por frame
# =============================================================================

def analizar_trayectorias(dcd_directory, ca_pdb, system, grupos=(1,2,3,4,5,12)):
    """
    Para cada trayectoria DCD en el directorio, calcula frame a frame:
        - Energía total y por grupo de fuerza
        - Q (parámetro de orden, grupo 12)
        - Rg (radio de giro)
        - RMSD respecto al frame inicial

    Guarda los resultados en archivos .txt para procesamiento con WHAM.
    """
    import os, glob
    dcd_files = glob.glob(os.path.join(dcd_directory, "*.dcd"))
    pdb_ref   = app.PDBFile(ca_pdb)

    for dcd_file in dcd_files:
        traj   = md.load_dcd(dcd_file, top=ca_pdb)
        rg     = md.compute_rg(traj)
        rmsd   = md.rmsd(traj, traj, frame=0)

        integrador = LangevinIntegrator(1 * kelvin, 0.5 / picosecond, 0.0005 * picoseconds)
        sim        = Simulation(pdb_ref.topology, system, integrador)
        sim.context.setPositions(pdb_ref.positions)

        base     = os.path.join(dcd_directory, os.path.basename(dcd_file))
        out_E    = open(base + '.txt',      'w')
        out_rg   = open(base + 'rg.txt',   'w')
        out_rmsd = open(base + 'rmsd.txt', 'w')

        for frame_idx in range(len(traj)):
            pos = traj[frame_idx].openmm_positions(0)
            sim.context.setPositions(pos)

            energias = []
            E_total  = 0.0
            E_Q      = 0.0
            for g in grupos:
                state = sim.context.getState(getEnergy=True, groups={g})
                E_g   = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
                energias.append(E_g)
                E_total += E_g
                if g == 12:
                    E_Q = E_g

            out_E.write(   f"{E_total} {E_Q} {E_Q}\n")
            out_rg.write(  f"{E_total} {E_Q} {rg[frame_idx]} {E_Q}\n")
            out_rmsd.write(f"{E_total} {E_Q} {rmsd[frame_idx]} {E_Q}\n")

        out_E.close(); out_rg.close(); out_rmsd.close()
        print(f"Analizado: {dcd_file}")


# =============================================================================
# SECCIÓN 13: VISUALIZACIÓN DE RESULTADOS
# =============================================================================

def graficar_energia_vs_tiempo(data_an="data.an"):
    """
    Grafica la energía potencial en función del paso de simulación.
    Útil para verificar convergencia y estabilidad.
    """
    df = pd.read_csv(data_an, sep=' ', header=None)
    # Columna 1 = paso, columna 2 = energía (formato: paso,energía,temperatura)
    plt.figure()
    plt.plot(df[0], df[1])
    plt.xlabel("Paso")
    plt.ylabel("Energía potencial (kJ/mol)")
    plt.title("Energía vs tiempo")
    plt.show()


def graficar_histogramas_Q(txt_directory):
    """
    Grafica histogramas de Q para cada ventana de umbrella.
    Las ventanas deben solaparse para que WHAM funcione correctamente.
    """
    import glob
    txt_files = glob.glob(os.path.join(txt_directory, "*.dcd.txt"))
    plt.figure()
    for f in txt_files:
        df = pd.read_csv(f, sep=r'\s+', header=None)
        df.columns = ['energia', 'Q', 'Qu']
        plt.hist(df['Q'], bins=50, alpha=0.5, label=os.path.basename(f))
    plt.xlabel("Q")
    plt.ylabel("Frecuencia")
    plt.title("Histogramas de Q por ventana de umbrella")
    plt.legend(fontsize=6)
    plt.show()
