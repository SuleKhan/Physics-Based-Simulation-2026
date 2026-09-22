#include "MassSpringApp.h"

#include <iostream>
#include <iomanip> 
#include <random>
#include <fstream>


static std::string vector3f_to_string(const Vector3F& v) {
    return std::to_string(v[0]) + " " +
           std::to_string(v[1]) + " " +
           std::to_string(v[2]);
}


class TestCase 
{
public:
    int m_method;
    double m_dt;
    double m_mass;
    double m_spring_stiffness;
    double m_spring_length;
    double m_damping;

    Vector3F m_gravity;
    Vector3F m_spring_start;

    Vector3F start_m_position;
    Vector3F start_m_velocity;
    double start_m_time;

    Vector3F expected_position;
    Vector3F expected_velocity;
    double expected_time;

    void save(std::ostream &out) const {

        out << std::fixed << std::setprecision(17);

        out << m_method << " "
            << m_dt << " "
            << m_mass << " "
            << m_spring_stiffness << " "
            << m_spring_length << " "
            << m_damping << " "
            << m_gravity[0] << " " << m_gravity[1] << " " << m_gravity[2] << " "
            << m_spring_start[0] << " " << m_spring_start[1] << " " << m_spring_start[2] << " "
            << start_m_position[0] << " " << start_m_position[1] << " " << start_m_position[2] << " "
            << start_m_velocity[0] << " " << start_m_velocity[1] << " " << start_m_velocity[2] << " "
            << start_m_time << " "
            << expected_position[0] << " " << expected_position[1] << " " << expected_position[2] << " "
            << expected_velocity[0] << " " << expected_velocity[1] << " " << expected_velocity[2] << " "
            << expected_time << std::endl;
    }

    bool load(const std::string& line) {
        std::istringstream in(line);
        in >> m_method
           >> m_dt
           >> m_mass
           >> m_spring_stiffness
           >> m_spring_length
           >> m_damping
           >> m_gravity[0] >> m_gravity[1] >> m_gravity[2]
           >> m_spring_start[0] >> m_spring_start[1] >> m_spring_start[2]
           >> start_m_position[0] >> start_m_position[1] >> start_m_position[2]
           >> start_m_velocity[0] >> start_m_velocity[1] >> start_m_velocity[2]
           >> start_m_time
           >> expected_position[0] >> expected_position[1] >> expected_position[2]
           >> expected_velocity[0] >> expected_velocity[1] >> expected_velocity[2]
           >> expected_time;

        return in.good();
    }

    void define_from_mass_spring_app(const MassSpringApp& app, 
                                     const Vector3F& start_m_vel,
                                     const Vector3F& start_m_pos,
                                     double start_m_t) {
        m_method = app.m_method;
        m_dt = app.m_dt;
        m_mass = app.m_mass;
        m_spring_stiffness = app.m_spring_stiffness;
        m_spring_length = app.m_spring_length;
        m_damping = app.m_damping;
        m_gravity = app.m_gravity;
        m_spring_start = app.m_spring_start;

        start_m_position = start_m_pos;
        start_m_velocity = start_m_vel;
        start_m_time = start_m_t;

        expected_position = app.m_position;
        expected_velocity = app.m_velocity;
        expected_time = app.m_time;
    }
};


std::string MassSpringApp::m_check() {

    std::string check_info;

    MassSpringApp original_state = *this;
    std::map<int, std::vector<TestCase>> tests_by_method;

    std::ifstream fin(PBS_DATA_DIR "/TestCases.dat");
    if (!fin) {
        return "ERROR: Could not open test file '" PBS_DATA_DIR "/TestCases.dat'!";
    }
    std::string line;
    std::vector<std::string> lines;
    while (std::getline(fin, line)) {
        lines.push_back(line);
    }
    fin.close();

    for (const auto& line : lines) {
        if (line.empty()) continue;
        TestCase test_case;
        if (!test_case.load(line + "\n")) {
            check_info += "Failed to load test line: " + line + "\n";
            continue;
        }
        tests_by_method[test_case.m_method].push_back(test_case);
    }
    

    for (const auto& pair : tests_by_method) {
        int method_id = pair.first;
        const std::vector<TestCase>& test_cases = pair.second;
        
        bool method_passed = true;
        
        for (const auto& test_case : test_cases) {
            m_damping = test_case.m_damping;
            m_dt = test_case.m_dt;
            m_gravity = test_case.m_gravity;
            m_mass = test_case.m_mass;
            m_method = test_case.m_method;
            m_spring_length = test_case.m_spring_length;
            m_spring_stiffness = test_case.m_spring_stiffness;
            m_spring_start = test_case.m_spring_start;
            m_time = test_case.start_m_time;
            m_position = test_case.start_m_position;
            m_velocity = test_case.start_m_velocity;

            step();

            bool correct = (m_position.isApprox(test_case.expected_position, 1e-6) &&
                            m_velocity.isApprox(test_case.expected_velocity, 1e-6));

            if (!correct) {
                check_info += "Method " + std::to_string(test_case.m_method) + ": FAILED\n";
                check_info += "  Expected next step position: " + vector3f_to_string(test_case.expected_position) + "\n";
                check_info += "  Actual next step position:   " + vector3f_to_string(m_position) + "\n";
                check_info += "  Position error: " + std::to_string((m_position - test_case.expected_position).norm()) + "\n";
                check_info += "  Expected next step velocity: " + vector3f_to_string(test_case.expected_velocity) + "\n";
                check_info += "  Actual next step velocity:   " + vector3f_to_string(m_velocity) + "\n";
                check_info += "  Velocity error: " + std::to_string((m_velocity - test_case.expected_velocity).norm()) + "\n";
                *this = original_state;
                return check_info;
            }
        }

        if (method_passed) {
            check_info += "Method " + std::to_string(method_id) + ": PASSED\n";
        }

    }

    *this = original_state;
    return check_info;
}


void MassSpringApp::createTests() {
    const std::string filename = "mass_spring_tests.txt";
    std::ofstream fout(filename, std::ios::app);
    if (!fout) {
        throw std::runtime_error("Could not open file!");
    }
    for (int method = 0; method <= 4; ++method) {
        m_method = method;
        m_position = Vector3F(0, -3.0, 0);
        m_velocity = Vector3F::Zero();
        for (int i = 0; i < 10; ++i) {
            step();
        }
        TestCase test_case;
        Vector3F start_m_pos = m_position;
        Vector3F start_m_vel = m_velocity;
        double start_m_t = m_time;
        step();
        test_case.define_from_mass_spring_app(*this, start_m_vel, start_m_pos, start_m_t);
        test_case.save(fout);
        std::cout << "Saved test case for method " << method << std::endl;
    }
}


void MassSpringApp::createRandomTests(int num_tests_per_method) {
    const std::string filename = "mass_spring_tests.txt";
    std::ofstream fout(filename, std::ios::app);
    if (!fout) {
        throw std::runtime_error("Could not open file!");
    }

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> pos_dist(-5.0, 5.0);
    std::uniform_real_distribution<double> vel_dist(-1.0, 1.0);
    std::uniform_real_distribution<double> time_dist(0.0, 10.0);

    for (int method = 0; method <= 4; ++method) {
        for (int i = 0; i < num_tests_per_method; ++i) {
            m_method = method;
            m_position = Vector3F(pos_dist(rng), pos_dist(rng), pos_dist(rng));
            m_velocity = Vector3F(vel_dist(rng), vel_dist(rng), vel_dist(rng));
            m_time = time_dist(rng);
            step();

            TestCase test_case;
            Vector3F start_m_pos = m_position;
            Vector3F start_m_vel = m_velocity;
            double start_m_t = m_time;
            step();
            test_case.define_from_mass_spring_app(*this, start_m_vel, start_m_pos, start_m_t);
            test_case.save(fout);
            std::cout << "Saved random test case for method " << method << std::endl;
        }
    }
}
