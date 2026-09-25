#ifndef SIMULATION_PARSER_H
#define SIMULATION_PARSER_H

#include <string>
#include <map>
#include <vector>

// ============================================================================
// Data Structures
// ============================================================================

/**
 * @brief Represents a simulation region with physical properties
 */
struct Region {
    double x1 = 0.0;          ///< Left boundary
    double x2 = 0.0;          ///< Right boundary
    double y1 = 0.0;          ///< Bottom boundary
    double y2 = 0.0;          ///< Top boundary
    double T = 0.0;           ///< Temperature
    double rho_0 = 0.0;       ///< Initial density
    double m = 0.0;           ///< Initial total mass
    bool do_output = false;   ///< Whether to output p, T data
};

/**
 * @brief Represents a wall boundary with constraint properties
 */
struct Wall {
    std::string type;         ///< Wall type (e.g., "fixed", "movable")
    double x1 = 0.0;          ///< Left boundary
    double x2 = 0.0;          ///< Right boundary
    double y1 = 0.0;          ///< Bottom boundary
    double y2 = 0.0;          ///< Top boundary
    bool lock_x = false;      ///< Lock x movement
    bool lock_y = false;      ///< Lock y movement
    bool lock_rot = false;    ///< Lock rotation
    double m = 0.0;           ///< Wall mass
    bool do_output = false;   ///< Whether to output position data
};

/**
 * @brief Represents a keystone point for piston definition
 */
struct Keystone {
    double x = 0.0;           ///< X coordinate
    double y = 0.0;           ///< Y coordinate
    double a = 0.0;           ///< Parameter a
    double t = 0.0;           ///< Parameter t
};

/**
 * @brief Represents a piston with keystones for motion definition
 */
struct Piston {
    double width = 0.0;       ///< Piston width
    double depth = 0.0;       ///< Piston depth
    double alpha = 0.0;       ///< Alpha parameter
    std::vector<Keystone> keystones;  ///< Keystone points
    bool loop = false;        ///< Whether to loop the motion
    bool do_output = false;   ///< Whether to output position data
};

/**
 * @brief Global simulation parameters
 */
struct SimulationParams {
    double sim_t = 0.0;       ///< Total simulation time
    double dt = 0.0;          ///< Base timestep
    int N_particle = 0;       ///< Number of total particles
};

// ============================================================================
// Parser Class Interface
// ============================================================================

/**
 * @class SimulationParser
 * @brief Parses simulation input files into structured data
 * 
 * This class reads simulation configuration files and stores the data
 * in accessible maps and structures for use in simulation code.
 * 
 * Example usage:
 * @code
 * SimulationParser parser;
 * if (parser.parse("simulation.input")) {
 *     double sim_time = parser.params.sim_t;
 *     auto& region = parser.regions["hot_zone"];
 *     // ... use parsed data
 * }
 * @endcode
 */
class SimulationParser {
public:
    // =========================================================================
    // Public Data Members (parsed results)
    // =========================================================================
    
    SimulationParams params;                          ///< Global simulation parameters
    std::map<std::string, Region> regions;            ///< Named regions
    std::map<std::string, Wall> walls;                ///< Named walls
    std::map<std::string, Piston> pistons;            ///< Named pistons

    // =========================================================================
    // Constructor/Destructor
    // =========================================================================
    
    SimulationParser() = default;
    ~SimulationParser() = default;

    // =========================================================================
    // Public Methods
    // =========================================================================
    
    /**
     * @brief Parse a simulation input file
     * @param filename Path to the input file
     * @return true if parsing succeeded, false otherwise
     */
    bool parse(const std::string& filename);

    /**
     * @brief Parse from a string directly (for testing/embedded configs)
     * @param content The input file content as a string
     * @return true if parsing succeeded, false otherwise
     */
    bool parseFromString(const std::string& content);

    /**
     * @brief Clear all parsed data
     */
    void clear();

    /**
     * @brief Check if a region exists
     * @param name Region name
     * @return true if region exists
     */
    bool hasRegion(const std::string& name) const;

    /**
     * @brief Check if a wall exists
     * @param name Wall name
     * @return true if wall exists
     */
    bool hasWall(const std::string& name) const;

    /**
     * @brief Check if a piston exists
     * @param name Piston name
     * @return true if piston exists
     */
    bool hasPiston(const std::string& name) const;

    /**
     * @brief Get the last error message
     * @return Error message string
     */
    std::string getLastError() const;

    /**
     * @brief Validate parsed data for consistency
     * @return true if data is valid, false otherwise
     */
    bool validate() const;

private:
    // =========================================================================
    // Private Methods (implementation details)
    // =========================================================================
    
    void trim(std::string& s);
    std::vector<std::string> tokenize(const std::string& line);
    bool toBool(const std::string& s);
    
    void parseGlobal(const std::string& line);
    void parseRegion(const std::string& name, const std::string& line);
    void parseWall(const std::string& name, const std::string& line);
    void parsePiston(const std::string& name, const std::string& line);
    void parseKeystone(const std::string& piston_name, const std::string& line);
    
    bool parseStream(std::istream& input);
    
    std::string last_error;
};

#endif // SIMULATION_PARSER_H
