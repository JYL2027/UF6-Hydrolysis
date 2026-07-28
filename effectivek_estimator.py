import math
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s"
)

logger = logging.getLogger(__name__)


# This piece of code/software calculates various approximations to the
# effective rate constant of UF6 hydrolysis based on user inputs.


def HK_pressure_decay_model(
    radius: float,
    nd: float,
    Erxn: float,
    Edes: float,
    Temp: float,
    nu: float,
    alpha: float
) -> float:
    """
    This function computes the effective rate constant for a Hertz-Knudsen
    pressure decay model of the form dP/dt = - k_eff PUF6.

    Parameters:
        radius: Radius of a water microdroplet in meters
        nd: Droplet number density in droplets/m^3
        Erxn: Reaction energy in kJ/mol
        Edes: Desorption energy in kJ/mol
        Temp: Temperature of the system/reaction in Kelvin
        nu: The attempt frequency in s^-1
        alpha: Accommodation coefficient, 0 to 1

    Output:
        keff: The effective rate constant in 1/s
    """

    logger.debug("Entering HK_pressure_decay_model()")
    logger.debug(
        "Inputs: radius=%s, nd=%s, Erxn=%s, Edes=%s, Temp=%s, nu=%s, alpha=%s",
        radius, nd, Erxn, Edes, Temp, nu, alpha
    )

    UF6mass = 5.845e-25

    # Gas constant in kJ/mol/K
    R = 0.008314

    # Boltzmann constant in J/K
    Boltz = 1.380649e-23

    try:
        krxn = nu * math.exp(-Erxn / (R * Temp))
        kdes = nu * math.exp(-Edes / (R * Temp))

        logger.debug("Computed krxn=%s", krxn)
        logger.debug("Computed kdes=%s", kdes)

    except Exception as e:
        logger.exception("Failed while computing krxn or kdes.")
        raise ValueError(f"Failed to compute krxn or kdes: {e}") from e

    try:
        term1 = (4 * math.pi * radius**2) * nd
        term2 = math.sqrt((Boltz * Temp) / (2 * math.pi * UF6mass)) * alpha
        term3 = krxn / (krxn + kdes)

        logger.debug("Computed term1=%s", term1)
        logger.debug("Computed term2=%s", term2)
        logger.debug("Computed term3=%s", term3)

        keff = term1 * term2 * term3

        logger.info("Computed HK pressure decay keff=%s 1/s", keff)

    except Exception as e:
        logger.exception("Failed while computing keff in HK_pressure_decay_model().")
        raise ValueError(f"Failed to compute keff: {e}") from e

    logger.debug("Exiting HK_pressure_decay_model()")

    return keff


def MESS_HK_preliminary_keff(
    keq: float,
    k2: float,
    radius: float,
    alpha: float,
    Temp: float,
    nu: float,
    Edes: float,
    CH2O: float
) -> tuple[float, float]:
    """
    This function computes the k_eff of the MESS preliminary model in two forms.

    One form is represented in the traditional sense of kinetics.
    The other is in a form representative of the pressure decay model.

    This model uses the Hertz-Knudsen approximation to the adsorption rate
    constant and the Arrhenius form of the desorption kinetics for estimating
    the desorption rate constant.

    rate = keff [H2O] [UF6]

    Parameters:
        keq: Equilibrium rate constant for dissociation calculated from MESS outputs
        k2: Forward reaction rate constant calculated from MESS outputs
        radius: Radius of a water microdroplet in meters
        alpha: Accommodation coefficient, 0 to 1
        Temp: Temperature of the system/reaction in Kelvin
        nu: Attempt frequency in s^-1
        Edes: Desorption energy in kJ/mol
        CH2O: Concentration of water in molecules/cm^3

    Output:
        keff_1: Effective rate constant of the form rate = keff [H2O] [UF6]
        keff_2: Effective rate constant of the form rate = keff [UF6]
    """

    logger.debug("Entering MESS_HK_preliminary_keff()")
    logger.debug(
        "Inputs: keq=%s, k2=%s, radius=%s, alpha=%s, Temp=%s, nu=%s, Edes=%s, CH2O=%s",
        keq, k2, radius, alpha, Temp, nu, Edes, CH2O
    )

    UF6mass = 5.845e-25

    # Gas constant in kJ/mol/K
    R = 0.008314

    # Boltzmann constant in J/K
    Boltz = 1.380649e-23

    try:
        kdes = nu * math.exp(-Edes / (R * Temp))

        kads = (
            (3 * alpha) / radius
        ) * math.sqrt((Boltz * Temp) / (2 * math.pi * UF6mass))

        logger.debug("Computed kdes=%s", kdes)
        logger.debug("Computed kads=%s", kads)

    except Exception as e:
        logger.exception("Failed while computing kdes or kads.")
        raise ValueError(f"Failed to compute kdes or kads: {e}") from e

    try:
        keff_1 = keq * k2 * (kads / kdes)
        keff_2 = keq * k2 * (kads / kdes) * CH2O

        logger.info("Computed MESS_HK keff_1=%s", keff_1)
        logger.info("Computed MESS_HK keff_2=%s", keff_2)

    except Exception as e:
        logger.exception("Failed while computing MESS_HK effective rate constants.")
        raise ValueError(f"Failed to compute MESS_HK effective rate constants: {e}") from e

    logger.debug("Exiting MESS_HK_preliminary_keff()")

    return keff_1, keff_2


def BEMM_method_flow_rate(p_inf: float, T_1: float, T_inf: float, velocity: float) -> float:
    '''
    This function solves for the evaporation flow rate using the nonlinearized BEMM.

    Parameters:
        p_inf: Pressure of the vapor (Pa)
        T_1: Temperature of the surface of the water micro-droplet (K)
        T_inf: Temperature of the vapor (K)
        velocity: Gas velocity outside the Knudsen layer
    
    Output:
        BEMM_j : The mass flux of the nonlinearized BEMM (kg/(m^2 s))
    
    '''
    logger.debug("Entering BEMM_method_flow_rate()")
    logger.debug(
        "Inputs: p_inf=%s Pa, T_1=%s K, T_inf=%s K, velocity=%s m/s",
        p_inf, T_1, T_inf, velocity
    )

    # molar mass of UF6s
    M_UF6 = 0.352
    Gas_constant = 8.314

    R = Gas_constant/M_UF6

    try:
        speed_ratio = velocity / math.sqrt(2 * R * T_inf)
        term1 = p_inf / math.sqrt(R * T_1)
        term2 = math.sqrt(2 / (T_inf / T_1))
        BEMM_j = term1 * term2 * speed_ratio
    except Exception as e:
        raise ValueError(f"{e}")

    logger.info("Computed BEMM mass flux=%s", BEMM_j)

    logger.debug("Exiting BEMM_method_flow_rate()")

    return BEMM_j


def BEMM_keff_from_flux(
    BEMM_j: float,
    p_inf: float,
    T_inf: float,
    radius: float,
    nd: float,
    reaction_probability: float
) -> float:
    '''
    This function converts the BEMM interfacial mass flux to an effective pressure-decay rate constant.

    Parameters:
        BEMM_j: BEMM mass flux.
            Units: kg/(m^2 s)

        p_inf: UF6 partial pressure in the gas phase.
            Units: Pa

        T_inf: Vapor/gas temperature outside the Knudsen layer.
            Units: K

        radius: Droplet radius.
            Units: m

        nd: Droplet number density.
            Units: droplets/m^3

        reaction_probability:
            Probabilty of reaction before desorbing

    Returns:
        keff: Effective pressure-decay rate constant.
            Units: 1/s
    '''

    logger.debug("Entering BEMM_keff_from_flux()")
    logger.debug(
        "Inputs: BEMM_j=%s, p_inf=%s, T_inf=%s, radius=%s, nd=%s, reaction_probability=%s",
        BEMM_j, p_inf, T_inf, radius, nd, reaction_probability
    )

    M_UF6 = 0.352  # kg/mol
    Gas_constant = 8.314462618  # J/(mol K)

    R = Gas_constant / M_UF6  # J/(kg K)

    try:
        area_per_volume = 4.0 * math.pi * radius**2 * nd

        keff = (
            area_per_volume
            * R
            * T_inf
            * BEMM_j
            / p_inf
            * reaction_probability
        )

    except Exception as e:
        logger.exception("Failed inside BEMM_keff_from_flux().")
        raise ValueError(f"{e}") from e

    logger.debug("Computed area_per_volume=%s 1/m", area_per_volume)
    logger.debug("Computed R=%s J/(kg K)", R)
    logger.info("Computed BEMM keff=%s 1/s", keff)

    logger.debug("Exiting BEMM_keff_from_flux()")

    return keff

def linearizedBEMM(alpha: float, T_1: float, radius: float, T_inf: float, p_inf: float, p_1: float, nd: float, reaction_probability: float ) -> float:
    '''
    This function returns the effective rate constant for the linearized BEMM. This is not as accurate as the normal BEMM and is a approximation that over estimates values. 

    Parameters: 
        alpha: accomodation coefficeint
        T_1: Temperature (K) of the water micro droplet surface 
        radius: Radius of the micro droplet (m)
        T_inf: Temperature (K) of the UF6 gas/vapor
        p_inf: Pressure of the Uf6 gas/vapor (Pa)
        p_1: Effective UF6 partial pressure at the water droplet surface, Pa
        nd: The droplet number density of the water micro droplet (droplets/m^3)
        reaction_probability: Probabilty of reaction before desorbing

    Returns:
        LinearBEMMK: Linearized BEMM effective rate constant 
    '''
    logger.debug("Entering linearizedBEMM()")
    logger.debug(
    "Inputs: alpha=%s, T_1=%s, radius=%s, T_inf=%s, p_inf=%s, p_1=%s, nd=%s, reaction_probability=%s",
    alpha, T_1, radius, T_inf, p_inf, p_1, nd, reaction_probability
    )

    M_UF6 = 0.352  # kg/mol
    Gas_constant = 8.314462618  # J/(mol K)
    
    R = Gas_constant / M_UF6  # J/(kg K)
    gamma = 1.67
    evaporation_flux = None

    # Calculating the evaporation flux
    try:
        if not 0.0 <= alpha <= 1.0:
            raise ValueError("alpha must be between 0 and 1.")

        if T_1 <= 0.0:
            raise ValueError("T_1 must be positive.")

        if T_inf <= 0.0:
            raise ValueError("T_inf must be positive.")

        if radius <= 0.0:
            raise ValueError("radius must be positive.")

        if p_inf <= 0.0:
            raise ValueError("p_inf must be positive.")

        if p_1 < 0.0:
            raise ValueError("p_1 must be nonnegative.")

        if nd < 0.0:
            raise ValueError("nd must be nonnegative.")

        if not 0.0 <= reaction_probability <= 1.0:
            raise ValueError("reaction_probability must be between 0 and 1.")

        # reverse the pressure as adsorption not evaporation
        delta_p = p_inf - p_1
        gamma_ratio = (gamma - 1) / gamma 
        term1 = alpha / (1 - (gamma_ratio * alpha))

        evaporation_flux = term1 * (delta_p / math.sqrt(2 * math.pi * R * T_1))
    except Exception as e:
        logger.exception("Failed inside linearizedBEMM().")
        raise ValueError(f"{e}") from e

    # Convert the flux to a first order rate constant 
    try:
        area_per_volume = 4.0 * math.pi * radius**2 * nd
    
        LinearBEMMK = evaporation_flux * area_per_volume* R* T_inf / p_inf* reaction_probability
                        
    except Exception as e: 
        raise ValueError(f"{e}") from e

    return LinearBEMMK

def solve(method: str, data: list) -> list:
    """
    This function takes user inputs and returns the requested effective
    rate constant or constants.

    Parameters:
        method: str
            Calculation method. Supported methods are:
            - 'MESS_HK'
            - 'HK_P'
            - 'BEMM'
            - 'LinearBEMM'

        data: list
            Input parameters corresponding to the selected method.

    Returns:
        Effective rate constant or constants as a list.
    """

    logger.debug("Entering solve()")
    logger.debug("Requested method=%s", method)
    logger.debug("Input data=%s", data)

    if method == "MESS_HK":
        logger.info("Solving using MESS_HK method.")

        try:
            keq, k2, radius, alpha, Temp, nu, Edes, CH2O = data
        except ValueError as e:
            logger.exception("Invalid data length for MESS_HK.")
            raise ValueError(
                "MESS_HK requires data = "
                "[keq, k2, radius, alpha, Temp, nu, Edes, CH2O]"
            ) from e

        logger.debug(
            "Unpacked MESS_HK data: keq=%s, k2=%s, radius=%s, alpha=%s, "
            "Temp=%s, nu=%s, Edes=%s, CH2O=%s",
            keq, k2, radius, alpha, Temp, nu, Edes, CH2O
        )

        if radius <= 0:
            logger.error("Invalid radius=%s. Radius must be positive.", radius)
            raise ValueError("radius must be positive.")

        if Temp <= 0:
            logger.error("Invalid Temp=%s. Temp must be positive Kelvin.", Temp)
            raise ValueError("Temp must be positive Kelvin.")

        if nu <= 0:
            logger.error("Invalid nu=%s. nu must be positive.", nu)
            raise ValueError("nu must be positive.")

        if not 0 <= alpha <= 1:
            logger.error("Invalid alpha=%s. Alpha must be between 0 and 1.", alpha)
            raise ValueError("alpha must be between 0 and 1.")

        if CH2O < 0:
            logger.error("Invalid CH2O=%s. CH2O must be nonnegative.", CH2O)
            raise ValueError("CH2O must be nonnegative.")

        keff1, keff2 = MESS_HK_preliminary_keff(
            keq, k2, radius, alpha, Temp, nu, Edes, CH2O
        )

        result = [keff1, keff2]

        logger.info("MESS_HK solve result=%s", result)
        logger.debug("Exiting solve()")

        return result

    elif method == "HK_P":
        logger.info("Solving using HK_P method.")

        try:
            radius, nd, Erxn, Edes, Temp, nu, alpha = data
        except ValueError as e:
            logger.exception("Invalid data length for HK_P.")
            raise ValueError(
                "HK_P requires data = "
                "[radius, nd, Erxn, Edes, Temp, nu, alpha]"
            ) from e

        logger.debug(
            "Unpacked HK_P data: radius=%s, nd=%s, Erxn=%s, Edes=%s, "
            "Temp=%s, nu=%s, alpha=%s",
            radius, nd, Erxn, Edes, Temp, nu, alpha
        )

        if radius <= 0:
            logger.error("Invalid radius=%s. Radius must be positive.", radius)
            raise ValueError("radius must be positive.")

        if nd < 0:
            logger.error("Invalid nd=%s. Droplet number density must be nonnegative.", nd)
            raise ValueError("droplet number density must be nonnegative.")

        if Temp <= 0:
            logger.error("Invalid Temp=%s. Temp must be positive Kelvin.", Temp)
            raise ValueError("Temp must be positive Kelvin.")

        if nu <= 0:
            logger.error("Invalid nu=%s. nu must be positive.", nu)
            raise ValueError("nu must be positive.")

        if not 0 <= alpha <= 1:
            logger.error("Invalid alpha=%s. Alpha must be between 0 and 1.", alpha)
            raise ValueError("alpha must be between 0 and 1.")

        keff = HK_pressure_decay_model(radius, nd, Erxn, Edes, Temp, nu, alpha)

        result = [keff]

        logger.info("HK_P solve result=%s", result)
        logger.debug("Exiting solve()")

        return result
    
    elif method == "BEMM":
        logger.info("Solving using BEMM method.")

        try:
            p_inf, T_1, T_inf, velocity, radius, nd, reaction_probability = data

        except (ValueError, TypeError) as e:
            logger.exception("Invalid data length for BEMM.")
            raise ValueError(
                "BEMM requires data = "
                "[p_inf, T_1, T_inf, velocity, radius, nd, reaction_probability]"
            ) from e

        logger.debug(
            "Unpacked BEMM data: p_inf=%s, T_1=%s, T_inf=%s, velocity=%s, "
            "radius=%s, nd=%s, reaction_probability=%s",
            p_inf, T_1, T_inf, velocity, radius, nd, reaction_probability
        )

        if radius <= 0:
            logger.error("Invalid radius=%s. Radius must be positive.", radius)
            raise ValueError("radius must be positive.")
        
        if nd < 0:
            logger.error("Invalid nd=%s. Droplet number density must be nonnegative.", nd)
            raise ValueError("droplet number density must be nonnegative.")
        
        if T_1 <= 0:
            logger.error("Invalid Temp=%s. Temp must be positive Kelvin at the droplet surface.", T_1)
            raise ValueError("Temp must be positive Kelvin.")

        if T_inf <= 0:
            logger.error("Invalid Temp=%s. Temp must be positive Kelvin at the vapor region.", T_inf)
            raise ValueError("Temp must be positive Kelvin.")
        
        if not 0 <= reaction_probability <= 1:
            logger.error("Invalid reaction_probability=%s. Reaction probability must be between 0 and 1.", reaction_probability)
            raise ValueError("reaction_probability must be between 0 and 1.")
        
        BEMM_flux = BEMM_method_flow_rate(
            p_inf=p_inf,
            T_1=T_1,
            T_inf=T_inf,
            velocity=velocity
        )

        logger.debug("Computed BEMM_flux=%s kg/(m^2 s)", BEMM_flux)

        result = BEMM_keff_from_flux(
            BEMM_j=BEMM_flux,
            p_inf=p_inf,
            T_inf=T_inf,
            radius=radius,
            nd=nd,
            reaction_probability=reaction_probability
        )

        logger.info("BEMM effective rate constant=%s 1/s", result)

        return [result]

    
    elif method == 'LinearBEMM':
        logger.info("Solving using Linear BEMM method.")
        
        try:
            alpha, T_1, radius, T_inf, p_inf, p_1, nd, reaction_probability = data
        
        except (ValueError, TypeError) as e:
            logger.exception("Invalid data length for Linear BEMM.")
            raise ValueError(
                    "Linear BEMM requires data = "
                "[alpha, T_1, radius, T_inf, p_inf, p_1, nd, reaction_probability]"
            ) from e
        
        logger.debug(
            "Unpacked Linear BEMM data: alpha=%s, T_1=%s, radius=%s, T_inf=%s, p_inf=%s, p_1=%s, "
            "nd=%s, reaction_probability=%s",
            alpha, T_1, radius, T_inf, p_inf, p_1, nd, reaction_probability
        )

        result = linearizedBEMM(alpha, T_1, radius, T_inf, p_inf, p_1, nd, reaction_probability)

        return [result]
    else:
        logger.error("Unknown method provided: %s", method)
        raise ValueError(f"Unknown method provided: {method}")

def main():
    logger.info("Starting UF6 hydrolysis effective rate constant calculator.")
    print("Available methods:")
    print("MESS_HK")
    print("HK_P")
    print("BEMM")
    print("LinearBEMM")

    method = input("What method would you like to use?\n").strip()

    data_input = input("Please input all the data necessary to compute the effective rate constant. Please enter comma-separated numeric values.\n" )

    try:
        data = [float(value.strip()) for value in data_input.split(",")]
        logger.debug("Parsed data input: %s", data)

    except ValueError as e:
        logger.exception("Failed to parse user input data.")
        raise ValueError(
            "All data values must be numeric and comma-separated."
        ) from e

    results = solve(method, data)

    print(f"The effective rate constant or constants with method {method} is/are {results} 1/s")

    logger.info("Finished.")



if __name__ == "__main__":
    main()
