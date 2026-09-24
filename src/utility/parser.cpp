#include "parser.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>

bool SimulationParser::parse(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Cannot open file " << filename << std::endl;
    return false;
  }
  
  std::string line;
  std::string current_section;
  std::string current_name;
  bool in_keystones = false;
  
  while (std::getline(file, line)) {
    trim(line);
    if (line.empty() || line[0] == '#') continue;
    
    if (line == "Region" || line == "Wall" || line == "Piston") {
      current_section = line;
      in_keystones = false;
      continue;
    }

    if (line == "{") {
      std::getline(file, line);
      trim(line);
      current_name = line;
      continue;
    }

    if (line == "keystones" && current_section == "Piston") {
      in_keystones = true;
      continue;
    }

    if (line == "}") {
      in_keystones = false;
      continue;
    }

    if (in_keystones && current_section == "Piston") {
      parseKeystone(current_name, line);
      continue;
    }

    if (current_section == "Region") {
      parseRegion(current_name, line);
    } else if (current_section == "Wall") {
      parseWall(current_name, line);
    } else if (current_section == "Piston") {
      parsePiston(current_name, line);
    } else {
      parseGlobal(line);
    }
  }

  return true;
}


void SimulationParser::trim(std::string& s) {
  s.erase(0, s.find_first_not_of(" \t\r\n"));
  s.erase(s.find_last_not_of(" \t\r\n") + 1);
}

std::vector<std::string> SimulationParser::tokenize(const std::string& line) {
  std::vector<std::string> tokens;
  std::istringstream iss(line);
  std::string token;
  while (iss >> token) {
    tokens.push_back(token);
  }
  return tokens;
}

bool SimulationParser::toBool(const std::string& s) {
  return (s == "T" || s == "true" || s == "1");
}

void SimulationParser::parseGlobal(const std::string& line) {
  auto tokens = tokenize(line);
  if (tokens.size() < 2) return;

  std::string key = tokens[0];
  std::string value = tokens[1];

  if (key == "sim_t") params.sim_t = std::stod(value);
  else if (key == "dt") params.dt = std::stod(value);
  else if (key == "N_particle") params.N_particle = std::stoi(value);
}

void SimulationParser::parseRegion(const std::string& name, const std::string& line) {
  auto tokens = tokenize(line);
  if (tokens.size() < 2) return;

  Region& r = regions[name];
  std::string key = tokens[0];
  std::string value = tokens[1];

  if (key == "x1") r.x1 = std::stod(value);
  else if (key == "x2") r.x2 = std::stod(value);
  else if (key == "y1") r.y1 = std::stod(value);
  else if (key == "y2") r.y2 = std::stod(value);
  else if (key == "T") r.T = std::stod(value);
  else if (key == "rho_0") r.rho_0 = std::stod(value);
  else if (key == "m") r.m = std::stod(value);
  else if (key == "do_output") r.do_output = toBool(value);
}

void SimulationParser::parseWall(const std::string& name, const std::string& line) {
  auto tokens = tokenize(line);
  if (tokens.size() < 2) return;

  Wall& w = walls[name];
  std::string key = tokens[0];
  std::string value = tokens[1];

  if (key == "type") w.type = value;
  else if (key == "x1") w.x1 = std::stod(value);
  else if (key == "x2") w.x2 = std::stod(value);
  else if (key == "y1") w.y1 = std::stod(value);
  else if (key == "y2") w.y2 = std::stod(value);
  else if (key == "lock_x") w.lock_x = toBool(value);
  else if (key == "lock_y") w.lock_y = toBool(value);
  else if (key == "lock_rot") w.lock_rot = toBool(value);
  else if (key == "m") w.m = std::stod(value);
  else if (key == "do_output") w.do_output = toBool(value);
}

void SimulationParser::parsePiston(const std::string& name, const std::string& line) {
  auto tokens = tokenize(line);
  if (tokens.size() < 2) return;

  Piston& p = pistons[name];
  std::string key = tokens[0];
  std::string value = tokens[1];

  if (key == "width") p.width = std::stod(value);
  else if (key == "depth") p.depth = std::stod(value);
  else if (key == "alpha") p.alpha = std::stod(value);
  else if (key == "loop") p.loop = toBool(value);
  else if (key == "do_output") p.do_output = toBool(value);
}

void SimulationParser::parseKeystone(const std::string& piston_name, const std::string& line) {
  auto tokens = tokenize(line);
  if (tokens.size() < 4) return;

  Keystone k;
  k.x = std::stod(tokens[0]);
  k.y = std::stod(tokens[1]);
  k.a = std::stod(tokens[2]);
  k.t = std::stod(tokens[3]);

  pistons[piston_name].keystones.push_back(k);
}

// ============================================================================
// Example Usage
// ============================================================================

int main(int argc, char* argv[]) {
    std::string filename = "simulation.input";
    if (argc > 1) filename = argv[1];

    SimulationParser parser;
    if (!parser.parse(filename)) {
        std::cerr << "Failed to parse input file." << std::endl;
        return 1;
    }

    // Output parsed data
    std::cout << "=== Simulation Parameters ===" << std::endl;
    std::cout << "sim_t: " << parser.params.sim_t << std::endl;
    std::cout << "dt: " << parser.params.dt << std::endl;
    std::cout << "N_particle: " << parser.params.N_particle << std::endl;

    std::cout << "\n=== Regions ===" << std::endl;
    for (const auto& [name, region] : parser.regions) {
        std::cout << "Region: " << name << std::endl;
        std::cout << "  Bounds: (" << region.x1 << ", " << region.x2 
                  << ") x (" << region.y1 << ", " << region.y2 << ")" << std::endl;
        std::cout << "  T: " << region.T << ", rho_0: " << region.rho_0 
                  << ", m: " << region.m << std::endl;
        std::cout << "  do_output: " << (region.do_output ? "true" : "false") << std::endl;
    }

    std::cout << "\n=== Walls ===" << std::endl;
    for (const auto& [name, wall] : parser.walls) {
        std::cout << "Wall: " << name << " (type: " << wall.type << ")" << std::endl;
        std::cout << "  Bounds: (" << wall.x1 << ", " << wall.x2 
                  << ") x (" << wall.y1 << ", " << wall.y2 << ")" << std::endl;
        std::cout << "  lock: x=" << wall.lock_x << " y=" << wall.lock_y 
                  << " rot=" << wall.lock_rot << ", m: " << wall.m << std::endl;
        std::cout << "  do_output: " << (wall.do_output ? "true" : "false") << std::endl;
    }

    std::cout << "\n=== Pistons ===" << std::endl;
    for (const auto& [name, piston] : parser.pistons) {
        std::cout << "Piston: " << name << std::endl;
        std::cout << "  width: " << piston.width << ", depth: " << piston.depth 
                  << ", alpha: " << piston.alpha << std::endl;
        std::cout << "  loop: " << (piston.loop ? "true" : "false") 
                  << ", do_output: " << (piston.do_output ? "true" : "false") << std::endl;
        std::cout << "  Keystones (" << piston.keystones.size() << "):" << std::endl;
        for (size_t i = 0; i < piston.keystones.size(); ++i) {
            const auto& k = piston.keystones[i];
            std::cout << "    [" << i << "] x=" << k.x << ", y=" << k.y 
                      << ", a=" << k.a << ", t=" << k.t << std::endl;
        }
    }

    return 0;
}
