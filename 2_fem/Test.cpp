#include "FEMApp.h"

#include <iostream>
#include <iomanip> 
#include <random>
#include <fstream>
#include <cassert>
#include <numbers>


class TestCase {
public:
    std::string test_name;
    std::vector<Element> elements;
    VectorXF input_x;
    std::vector<Pin> pins;
    double eps;
    bool use_gravity;

    std::vector<Matrix2F> expected_F;
    std::vector<double> expected_energies;
    int expected_pin_index;
    VectorXF expected_gradient;
    TripletListF expected_hessians;
    MatrixXF expected_H;

    void save(std::ostream &out) const {
        out << std::fixed << std::setprecision(17);
        out << test_name << " | ";
        out << elements.size() << " | ";
        for (const auto& elem : elements) {
            out << elem.vertices[0] << " " << elem.vertices[1] << " " << elem.vertices[2] << " "
                << elem.area << " "
                << elem.DX_inv(0,0) << " " << elem.DX_inv(0,1) << " "
                << elem.DX_inv(1,0) << " " << elem.DX_inv(1,1) << " | ";
        }
        
        if (test_name != "testInitElement") {
            out << input_x.size() << " ";
            for (int i = 0; i < input_x.size(); ++i) {
                out << input_x(i) << " ";
            }
            out << "| ";
        }

        if (test_name == "testDeformationGradient") {
            assert(elements.size() == expected_F.size());
            for (const auto& F : expected_F) {
                out << F(0,0) << " " << F(0,1) << " "
                    << F(1,0) << " " << F(1,1) << " | ";
            }
        }


        if (test_name == "testAddPinEnergy" ||
            test_name == "testAddPinGradient" ||
            test_name == "testAddPinHessian" || 
            test_name == "testComputeGradientFD" ||
            test_name == "testComputeHessianFD" || 
            test_name == "testAddSystemEnergy" ||
            test_name == "testAddSystemGradient" || 
            test_name == "testAddSystemHessian") {
            out << pins.size() << " ";
            for (const auto& pin : pins) {
                out << pin.node << " "
                    << pin.pos.x() << " " << pin.pos.y() << " "
                    << pin.k << " ";
            }
            out << "| ";
        }

        if (test_name == "testAddTriangleEnergy" || test_name == "testAddPinEnergy" || test_name == "testAddGravityEnergy" || test_name == "testAddSystemEnergy") {
            if (test_name == "testAddPinEnergy") {
                assert(pins.size() == expected_energies.size());
            } else if (test_name == "testAddTriangleEnergy") {
                assert(elements.size() == expected_energies.size());
            } else if (test_name == "testAddGravityEnergy" || test_name == "testAddSystemEnergy") {
                assert(1 == expected_energies.size());
            }

            for (const auto& energy : expected_energies) {
                out << energy << " ";
            }
            out << "| ";
        }

        if (test_name == "testAddPinGradient" || test_name == "testAddGravityGradient" || test_name == "testComputeGradientFD" || test_name == "testAddSystemGradient") {
            assert(input_x.size() == expected_gradient.size());
            if (test_name == "testAddPinGradient") {
                out << expected_pin_index << " ";
            }
            for (int i = 0; i < expected_gradient.size(); ++i) {
                out << expected_gradient(i) << " ";
            }
            out << "| ";
        }

        if (test_name == "testAddPinHessian" || test_name == "testAddGravityHessian" || test_name == "testAddSystemHessian") {
            if (test_name == "testAddPinHessian") {
                out << expected_pin_index << " ";
            }
            out << expected_hessians.size() << " ";
            for (const auto& triplet : expected_hessians) {
                out << triplet.col() << " " << triplet.row() << " " << triplet.value() << " ";
            }
            out << "| ";
        }

        if (test_name == "testComputeHessianFD" || test_name == "testComputeGradientFD") {
            out << eps << " ";
            out << "| ";
        }

        if (test_name == "testComputeHessianFD") {
            assert(expected_H.rows() == input_x.size() && expected_H.cols() == input_x.size());
            out << expected_H.rows() << " " << expected_H.cols() << " ";
            for (int r = 0; r < expected_H.rows(); ++r) {
                for (int c = 0; c < expected_H.cols(); ++c) {
                    out << expected_H(r,c) << " ";
                }
            }
            out << "| ";
        }

        if (test_name == "testAddSystemEnergy") {
            out << use_gravity << " ";
            out << "| ";
        }

        out << std::endl;
    }

    bool load(const std::string& line) {
        std::istringstream in(line);
        std::string buf;
        size_t num_elements;
        in >> test_name;
        in >> buf;
        in >> num_elements;
        in >> buf;
        elements.resize(num_elements);
        for (size_t i = 0; i < num_elements; ++i) {
            in >> elements[i].vertices[0] >> elements[i].vertices[1] >> elements[i].vertices[2]
               >> elements[i].area
               >> elements[i].DX_inv(0,0) >> elements[i].DX_inv(0,1)
               >> elements[i].DX_inv(1,0) >> elements[i].DX_inv(1,1);
            in >> buf;
        }

        if (test_name != "testInitElement") {  
            input_x.resize(num_elements);
            size_t x_size;
            in >> x_size;
            input_x.resize(x_size);
            for (size_t j = 0; j < x_size; ++j) {
                in >> input_x(j);
            }
            in >> buf;
        }

        if (test_name == "testDeformationGradient") {
            expected_F.resize(num_elements);
            for (size_t i = 0; i < num_elements; ++i) {
                in >> expected_F[i](0,0) >> expected_F[i](0,1)
                >> expected_F[i](1,0) >> expected_F[i](1,1);
                in >> buf;
            }
        }

        if (test_name == "testAddPinEnergy" ||
            test_name == "testAddPinGradient" ||
            test_name == "testAddPinHessian" || 
            test_name == "testComputeGradientFD" ||
            test_name == "testComputeHessianFD" || 
            test_name == "testAddSystemEnergy" ||
            test_name == "testAddSystemGradient" || 
            test_name == "testAddSystemHessian") {
            size_t num_pins;
            in >> num_pins;
            pins.resize(num_pins);
            for (size_t i = 0; i < num_pins; ++i) {
                in >> pins[i].node
                   >> pins[i].pos.x() >> pins[i].pos.y()
                   >> pins[i].k;
            }
            in >> buf;
        }

        if (test_name == "testAddTriangleEnergy" || test_name == "testAddPinEnergy" || test_name == "testAddGravityEnergy" || test_name == "testAddSystemEnergy") {
            if (test_name == "testAddPinEnergy") {
                expected_energies.resize(pins.size());
            } else if (test_name == "testAddTriangleEnergy") {
                expected_energies.resize(num_elements);
            } else if (test_name == "testAddGravityEnergy" || test_name == "testAddSystemEnergy") {
                expected_energies.resize(1);
            }
            for (size_t i = 0; i < expected_energies.size(); ++i) {
                in >> expected_energies[i];
            }
            in >> buf;
        }

        if (test_name == "testAddPinGradient" || test_name == "testAddGravityGradient" || test_name == "testComputeGradientFD" || test_name == "testAddSystemGradient") {
            if (test_name == "testAddPinGradient") {
                in >> expected_pin_index;
            }
            expected_gradient.resize(input_x.size());
            for (int i = 0; i < expected_gradient.size(); ++i) {
                in >> expected_gradient(i);
            }
            in >> buf;
        }

        if (test_name == "testAddPinHessian" || test_name == "testAddGravityHessian" || test_name == "testAddSystemHessian") {
            if (test_name == "testAddPinHessian") {
                in >> expected_pin_index;
            }
            size_t num_hessians;
            in >> num_hessians;
            expected_hessians.resize(num_hessians);
            for (size_t i = 0; i < expected_hessians.size(); ++i) {
                int col, row;
                double value;
                in >> col >> row >> value;
                expected_hessians[i] = TripletF(row, col, value);
            }
            in >> buf;
        }

        if (test_name == "testComputeHessianFD" || test_name == "testComputeGradientFD") {
            in >> eps;
            in >> buf;
        }

        if (test_name == "testComputeHessianFD") {
            int rows, cols;
            in >> rows >> cols;
            expected_H.resize(rows, cols);
            for (int r = 0; r < rows; ++r) {
                for (int c = 0; c < cols; ++c) {
                    in >> expected_H(r,c);
                }
            }
            in >> buf;
        }

        if (test_name == "testAddSystemEnergy") {
            in >> use_gravity;
            in >> buf;
        }

        return in.good();
    }

};


std::string testAddSystemEnergy(FEMApp &app, std::vector<TestCase> &test_cases) {
    for (auto &test_case : test_cases) {
        app.elements = test_case.elements;
        app.pins = test_case.pins;
        app.use_gravity = test_case.use_gravity;
        VectorXF &x = test_case.input_x;
        double E_test = 0.0;
        double E_expected = test_case.expected_energies[0];
        app.addSystemEnergy(x, E_test);
        if (std::abs(E_test - E_expected) > 1e-6) {
            std::string output = "testAddSystemEnergy FAILED: energy does not match\n";
            output += "  Computed energy: " + std::to_string(E_test) + ", Expected energy: " + std::to_string(E_expected) + "\n";
            return output; 
        } 
    }
    return "testAddSystemEnergy PASSED for all test cases.\n";
}


std::string testAddSystemGradient(FEMApp &app, std::vector<TestCase> &test_cases) {
    for (auto &test_case : test_cases) {
        app.elements = test_case.elements;
        app.use_gravity = test_case.use_gravity;
        app.pins = test_case.pins;
        VectorXF &x = test_case.input_x;
        VectorXF g_test = VectorXF::Zero(x.size());
        VectorXF g_expected = test_case.expected_gradient;
        app.addSystemGradient(x, g_test);
        if ((g_test - g_expected).norm() > 1e-6) {
            std::string output = "testAddSystemGradient FAILED: gradient does not match\n";
            output += "  Computed gradient: ";
            for (int j = 0; j < g_test.size(); ++j) {
                output += std::to_string(g_test(j)) + " ";
            }
            output += "\n";
            output += "  Expected gradient: ";
            for (int j = 0; j < g_expected.size(); ++j) {
                output += std::to_string(g_expected(j)) + " ";
            }
            output += "\n";
            return output; 
        } 
    }
    return "testAddSystemGradient PASSED for all test cases.\n";
}


std::string testAddSystemHessian(FEMApp &app, std::vector<TestCase> &test_cases) {
    for (auto &test_case : test_cases) {
        app.elements = test_case.elements;
        app.use_gravity = test_case.use_gravity;
        app.pins = test_case.pins;
        VectorXF &x = test_case.input_x;
        TripletListF hessian_triplets_test;
        TripletListF hessian_triplets_expected = test_case.expected_hessians;
        app.addSystemHessian(x, hessian_triplets_test);

        if (hessian_triplets_test.size() != hessian_triplets_expected.size()) {
            std::string output = "testAddSystemHessian FAILED: number of triplets does not match\n";
            output += "  Computed triplets: " + std::to_string(hessian_triplets_test.size()) + "\n";
            output += "  Expected triplets: " + std::to_string(hessian_triplets_expected.size()) + "\n";
            return output;
        }

        for (size_t i = 0; i < hessian_triplets_test.size(); ++i) {
            if (hessian_triplets_test[i].row() != hessian_triplets_expected[i].row() ||
                hessian_triplets_test[i].col() != hessian_triplets_expected[i].col() ||
                std::abs(hessian_triplets_test[i].value() - hessian_triplets_expected[i].value()) > 1e-6) {
                std::string output = "testAddSystemHessian FAILED: triplet " + std::to_string(i) + " does not match\n";
                output += "  Computed triplet: (" + std::to_string(hessian_triplets_test[i].row()) + ", "
                        + std::to_string(hessian_triplets_test[i].col()) + ", "
                        + std::to_string(hessian_triplets_test[i].value()) + ")\n";
                output += "  Expected triplet: (" + std::to_string(hessian_triplets_expected[i].row()) + ", "
                        + std::to_string(hessian_triplets_expected[i].col()) + ", "
                        + std::to_string(hessian_triplets_expected[i].value()) + ")\n";
                return output;
            }
        }
    }

    return "testAddSystemHessian PASSED for all test cases.\n";
}


std::string testComputeGradientFD(FEMApp &app, std::vector<TestCase> &test_cases) {
    for (auto &test_case : test_cases) {
        app.elements = test_case.elements;
        app.pins = test_case.pins;
        VectorXF &x = test_case.input_x;
        double eps = test_case.eps;
        VectorXF g_test = VectorXF::Zero(x.size());
        VectorXF g_expected = test_case.expected_gradient;
        app.computeGradientFD(x, g_test, eps);
        if ((g_test - g_expected).norm() > 1e-6) {
            std::string output = "testComputeGradientFD FAILED: gradient does not match\n";
            output += "  Computed gradient: ";
            for (int j = 0; j < g_test.size(); ++j) {
                output += std::to_string(g_test(j)) + " ";
            }
            output += "\n";
            output += "  Expected gradient: ";
            for (int j = 0; j < g_expected.size(); ++j) {
                output += std::to_string(g_expected(j)) + " ";
            }
            output += "\n";
            return output; 
        } 
    }
 
    return "testComputeGradientFD PASSED for all test cases.\n";
}


std::string testComputeHessianFD(FEMApp &app, std::vector<TestCase> &test_cases) {
    for (auto &test_case : test_cases) {
        app.elements = test_case.elements;
        app.pins = test_case.pins;
        VectorXF &x = test_case.input_x;
        double eps = test_case.eps;
        MatrixXF H_test = MatrixXF::Zero(x.size(), x.size());
        MatrixXF H_expected = test_case.expected_H;
        app.computeHessianFD(x, H_test, eps);
        if ((H_test - H_expected).norm() > 1e-6) {
            std::string output = "testComputeHessianFD FAILED: Hessian does not match\n";
            output += "  Difference:\n";
            for (int r = 0; r < H_test.rows(); ++r) {
                for (int c = 0; c < H_test.cols(); ++c) {
                    output += std::to_string(H_test(r,c) - H_expected(r,c)) + " ";
                }
                output += "\n";
            }
            return output;
        }
    }
 
    return "testComputeHessianFD PASSED for all test cases.\n";
}


std::string testAddGravityEnergy(FEMApp &app, std::vector<TestCase> &test_cases) {
    for (auto &test_case : test_cases) {
        app.elements = test_case.elements;
        VectorXF &x = test_case.input_x;
        double E_test = 0.0;
        double E_expected = test_case.expected_energies[0];
        app.addGravityEnergy(x, E_test);
        if (std::abs(E_test - E_expected) > 1e-6) {
            std::string output = "testAddGravityEnergy FAILED: energy does not match\n";
            output += "  Computed energy: " + std::to_string(E_test) + ", Expected energy: " + std::to_string(E_expected) + "\n";
            return output; 
        } 
    }
    return "testAddGravityEnergy PASSED for all test cases.\n";
}   


std::string testAddGravityGradient(FEMApp &app, std::vector<TestCase> &test_cases) {
    for (auto &test_case : test_cases) {
        app.elements = test_case.elements;
        VectorXF &x = test_case.input_x;
        VectorXF g_test = VectorXF::Zero(x.size());
        VectorXF g_expected = test_case.expected_gradient;
        app.addGravityGradient(x, g_test);
        if ((g_test - g_expected).norm() > 1e-6) {
            std::string output = "testAddGravityGradient FAILED: gradient does not match\n";
            output += "  Computed gradient: ";
            for (int j = 0; j < g_test.size(); ++j) {
                output += std::to_string(g_test(j)) + " ";
            }
            output += "\n";
            output += "  Expected gradient: ";
            for (int j = 0; j < g_expected.size(); ++j) {
                output += std::to_string(g_expected(j)) + " ";
            }
            output += "\n";
            return output; 
        } 
    }
 
    return "testAddGravityGradient PASSED for all test cases.\n";
}


std::string testAddGravityHessian(FEMApp &app, std::vector<TestCase> &test_cases) {
    for (auto &test_case : test_cases) {
        app.elements = test_case.elements;
        VectorXF &x = test_case.input_x;
        TripletListF hessian_triplets_test;
        TripletListF hessian_triplets_expected = test_case.expected_hessians;
        app.addGravityHessian(x, hessian_triplets_test);

        if (hessian_triplets_test.size() != hessian_triplets_expected.size()) {
            std::string output = "testAddGravityHessian FAILED: number of triplets does not match\n";
            output += "  Computed triplets: " + std::to_string(hessian_triplets_test.size()) + "\n";
            output += "  Expected triplets: " + std::to_string(hessian_triplets_expected.size()) + "\n";
            return output; 
        }

        for (size_t i = 0; i < hessian_triplets_test.size(); ++i) {
            const auto& triplet_test = hessian_triplets_test[i];
            const auto& triplet_expected = hessian_triplets_expected[i];
            if (triplet_test.row() != triplet_expected.row() ||
                triplet_test.col() != triplet_expected.col() ||
                std::abs(triplet_test.value() - triplet_expected.value()) > 1e-6) {
                std::string output = "testAddGravityHessian FAILED: triplet does not match at index " + std::to_string(i) + "\n";
                output += "  Computed triplet: (" + std::to_string(triplet_test.row()) + ", " + std::to_string(triplet_test.col()) + ", " + std::to_string(triplet_test.value()) + ")\n";
                output += "  Expected triplet: (" + std::to_string(triplet_expected.row()) + ", " + std::to_string(triplet_expected.col()) + ", " + std::to_string(triplet_expected.value()) + ")\n";
                return output; 
            }
        }
    }
    return "testAddGravityHessian PASSED for all test cases.\n";
}


std::string testAddPinEnergy(FEMApp &app, std::vector<TestCase> &test_cases) {
    for (auto &test_case : test_cases) {
        app.elements = test_case.elements;
        app.pins = test_case.pins;
        VectorXF &x = test_case.input_x;
        for (size_t i = 0; i < app.pins.size(); ++i) {
            double E_test = 0.0;
            double E_expected = test_case.expected_energies[i];
            app.addPinEnergy(x, i, E_test);
            if (std::abs(E_test - E_expected) > 1e-6) {
                std::string output = "testAddPinEnergy FAILED: energy does not match for pin " + std::to_string(i) + "\n";
                output += "  Pin node index: " + std::to_string(app.pins[i].node) + "\n";
                output += "  Computed energy: " + std::to_string(E_test) + ", Expected energy: " + std::to_string(E_expected) + "\n";
                return output; 
            } 
        }
    }
    return "testAddPinEnergy PASSED for all test cases.\n";
}


std::string testAddPinGradient(FEMApp &app, std::vector<TestCase> &test_cases) {
    for (auto &test_case : test_cases) {
        app.elements = test_case.elements;
        app.pins = test_case.pins;
        VectorXF &x = test_case.input_x;
        VectorXF g_test = VectorXF::Zero(x.size());
        VectorXF g_expected = test_case.expected_gradient;
        size_t pin_index = 0;
        app.addPinGradient(x, pin_index, g_test);
        if ((g_test - g_expected).norm() > 1e-6) {
            std::string output = "testAddPinGradient FAILED: gradient does not match\n";
            output += "  Computed gradient: ";
            for (int j = 0; j < g_test.size(); ++j) {
                output += std::to_string(g_test(j)) + " ";
            }
            output += "\n";
            output += "  Expected gradient: ";
            for (int j = 0; j < g_expected.size(); ++j) {
                output += std::to_string(g_expected(j)) + " ";
            }
            output += "\n";
            return output; 
        } 
    }
    return "testAddPinGradient PASSED for all test cases.\n";
}


std::string testAddPinHessian(FEMApp &app, std::vector<TestCase> &test_cases) {
    for (auto &test_case : test_cases) {
        app.elements = test_case.elements;
        app.pins = test_case.pins;
        VectorXF &x = test_case.input_x;
        TripletListF hessian_triplets_test;
        TripletListF hessian_triplets_expected = test_case.expected_hessians;
        size_t pin_index = 0;
        app.addPinHessian(x, pin_index, hessian_triplets_test);

        if (hessian_triplets_test.size() != hessian_triplets_expected.size()) {
            std::string output = "testAddPinHessian FAILED: number of triplets does not match\n";
            output += "  Computed triplets: " + std::to_string(hessian_triplets_test.size()) + "\n";
            output += "  Expected triplets: " + std::to_string(hessian_triplets_expected.size()) + "\n";
            return output; 
        }

        for (size_t i = 0; i < hessian_triplets_test.size(); ++i) {
            const auto& triplet_test = hessian_triplets_test[i];
            const auto& triplet_expected = hessian_triplets_expected[i];
            if (triplet_test.row() != triplet_expected.row() ||
                triplet_test.col() != triplet_expected.col() ||
                std::abs(triplet_test.value() - triplet_expected.value()) > 1e-6) {
                std::string output = "testAddPinHessian FAILED: triplet does not match at index " + std::to_string(i) + "\n";
                output += "  Computed triplet: (" + std::to_string(triplet_test.row()) + ", " + std::to_string(triplet_test.col()) + ", " + std::to_string(triplet_test.value()) + ")\n";
                output += "  Expected triplet: (" + std::to_string(triplet_expected.row()) + ", " + std::to_string(triplet_expected.col()) + ", " + std::to_string(triplet_expected.value()) + ")\n";
                return output; 
            }
        }
    }
    return "testAddPinHessian PASSED for all test cases.\n";
}


std::string testAddTriangleEnergy(FEMApp &app, std::vector<TestCase> &test_cases) {
    for (auto &test_case : test_cases) {
        app.elements = test_case.elements;
        VectorXF &x = test_case.input_x;
        for (size_t i = 0; i < app.elements.size(); ++i) {
            double E_test = 0.0;
            double E_expected = test_case.expected_energies[i];
            app.addTriangleEnergy(x, i, E_test);
            if (std::abs(E_test - E_expected) > 1e-6) {
                std::string output = "testAddTriangleEnergy FAILED: energy does not match for element " + std::to_string(i) + "\n";
                output += "  Vertex indices: " + std::to_string(app.elements[i].vertices[0]) + ", " + std::to_string(app.elements[i].vertices[1]) + ", " + std::to_string(app.elements[i].vertices[2]) + "\n";
                output += "  Computed energy: " + std::to_string(E_test) + ", Expected energy: " + std::to_string(E_expected) + "\n";
                return output; 
            } 
        }
    }
    return "testAddTriangleEnergy PASSED for all test cases.\n";
}


std::string testDeformationGradient(FEMApp &app, std::vector<TestCase> &test_cases) {
    for (auto &test_case : test_cases) {
        app.elements = test_case.elements;
        VectorXF &x = test_case.input_x;
        for (size_t i = 0; i < app.elements.size(); ++i) {
            Matrix2F F = app.deformationGradient(x, i);
            const Matrix2F &F_expected = test_case.expected_F[i];
            if ((F - F_expected).norm() > 1e-6) {
                std::string output = "testDeformationGradient FAILED: F does not match for element " + std::to_string(i) + "\n";
                output += "  Vertex indices: " + std::to_string(app.elements[i].vertices[0]) + ", " + std::to_string(app.elements[i].vertices[1]) + ", " + std::to_string(app.elements[i].vertices[2]) + "\n";
                output += "  Computed F: " + std::to_string(F(0,0)) + " " + std::to_string(F(0,1)) + " " + std::to_string(F(1,0)) + " " + std::to_string(F(1,1)) + "\n";
                output += "  Expected F: " + std::to_string(F_expected(0,0)) + " " + std::to_string(F_expected(0,1)) + " " + std::to_string(F_expected(1,0)) + " " + std::to_string(F_expected(1,1)) + "\n";
                return output;
            }
        }
    }

    return "testDeformationGradient PASSED for all test cases.\n";
}


std::string testInitElement(FEMApp &app, std::vector<TestCase> &test_cases) {

    for (auto &test_case : test_cases) {
        app.elements = test_case.elements;
        for (auto &elem : app.elements) {
            elem.area = 0.0;
            elem.DX_inv.setZero();
        }
        for (size_t i = 0; i < app.elements.size(); ++i) {
            app.initElement(i);
        }
        for (size_t i = 0; i < app.elements.size(); ++i) {
            const Element &elem1 = app.elements[i];
            const Element &elem2 = test_case.elements[i];
            if (std::abs(elem1.area - elem2.area) > 1e-6) {
                std::string output = "testInitElement FAILED: area does not match for element " + std::to_string(i) + "\n";
                output += "  Vertex indices: " + std::to_string(elem1.vertices[0]) + ", " + std::to_string(elem1.vertices[1]) + ", " + std::to_string(elem1.vertices[2]) + "\n";
                output += "  Computed area: " + std::to_string(elem1.area) + ", Expected area: " + std::to_string(elem2.area) + "\n";
                return output;

            }
            if ((elem1.DX_inv - elem2.DX_inv).norm() > 1e-6) {
                std::string output = "testInitElement FAILED: DX_inv does not match for element " + std::to_string(i) + "\n";
                output += "  Vertex indices: " + std::to_string(elem1.vertices[0]) + ", " + std::to_string(elem1.vertices[1]) + ", " + std::to_string(elem1.vertices[2]) + "\n";
                output += "  Computed DX_inv: " + std::to_string(elem1.DX_inv(0,0)) + " " + std::to_string(elem1.DX_inv(0,1)) + " " + std::to_string(elem1.DX_inv(1,0)) + " " + std::to_string(elem1.DX_inv(1,1)) + "\n";
                output += "  Expected DX_inv: " + std::to_string(elem2.DX_inv(0,0)) + " " + std::to_string(elem2.DX_inv(0,1)) + " " + std::to_string(elem2.DX_inv(1,0)) + " " + std::to_string(elem2.DX_inv(1,1)) + "\n";
                return output;
            }
        }
        std::string check_info;
    }

    return "testInitElement PASSED for all test cases.\n";
}


std::map<std::string, std::function<std::string(FEMApp&, std::vector<TestCase>&)>> tests = {
    {"testAddSystemEnergy", testAddSystemEnergy},
    {"testAddSystemGradient", testAddSystemGradient},
    {"testAddSystemHessian", testAddSystemHessian},
    {"testComputeGradientFD", testComputeGradientFD},
    {"testComputeHessianFD", testComputeHessianFD},
    {"testAddGravityEnergy", testAddGravityEnergy},
    {"testAddGravityGradient", testAddGravityGradient},
    {"testAddGravityHessian", testAddGravityHessian},
    {"testAddPinEnergy", testAddPinEnergy},
    {"testAddPinGradient", testAddPinGradient},
    {"testAddPinHessian", testAddPinHessian},
    {"testAddTriangleEnergy", testAddTriangleEnergy},
    {"testDeformationGradient", testDeformationGradient},
    {"testInitElement", testInitElement},
};


void generateFEMTests(FEMApp &original_app, bool use_original_app = false) {
    FEMApp app;
    if (use_original_app) {
        app = original_app;
    } else {
        app = FEMApp();
        app.nx = 5;
        app.ny = 5;
        app.resetSimulation();
    }
    
    std::vector<TestCase> test_cases;

    TestCase init_elem_test_case;
    init_elem_test_case.test_name = "testInitElement";
    init_elem_test_case.elements = app.elements;
    test_cases.push_back(init_elem_test_case);

    TestCase def_grad_test_case;
    def_grad_test_case.test_name = "testDeformationGradient";
    def_grad_test_case.elements = app.elements;
    
    TestCase add_tri_energy_test_case;
    add_tri_energy_test_case.test_name = "testAddTriangleEnergy";
    add_tri_energy_test_case.elements = app.elements;

    TestCase add_pin_energy_test_case;
    add_pin_energy_test_case.test_name = "testAddPinEnergy";
    add_pin_energy_test_case.elements = app.elements;

    TestCase add_pin_grad_test_case;
    add_pin_grad_test_case.test_name = "testAddPinGradient";
    add_pin_grad_test_case.elements = app.elements;

    TestCase add_pin_hess_test_case;
    add_pin_hess_test_case.test_name = "testAddPinHessian";
    add_pin_hess_test_case.elements = app.elements;

    TestCase add_grav_energy_test_case;
    add_grav_energy_test_case.test_name = "testAddGravityEnergy";
    add_grav_energy_test_case.elements = app.elements;
    
    TestCase add_grav_grad_test_case;
    add_grav_grad_test_case.test_name = "testAddGravityGradient";
    add_grav_grad_test_case.elements = app.elements;
    
    TestCase add_grav_hess_test_case;
    add_grav_hess_test_case.test_name = "testAddGravityHessian";
    add_grav_hess_test_case.elements = app.elements;

    TestCase compute_grad_fd_test_case;
    compute_grad_fd_test_case.test_name = "testComputeGradientFD";
    compute_grad_fd_test_case.elements = app.elements;

    TestCase compute_hess_fd_test_case;
    compute_hess_fd_test_case.test_name = "testComputeHessianFD";
    compute_hess_fd_test_case.elements = app.elements;

    TestCase add_sys_energy_test_case;
    add_sys_energy_test_case.test_name = "testAddSystemEnergy";
    add_sys_energy_test_case.elements = app.elements;

    TestCase add_sys_grad_test_case;
    add_sys_grad_test_case.test_name = "testAddSystemGradient";
    add_sys_grad_test_case.elements = app.elements;

    TestCase add_sys_hess_test_case;
    add_sys_hess_test_case.test_name = "testAddSystemHessian";
    add_sys_hess_test_case.elements = app.elements;

    if (!use_original_app) {
        app.use_dynamics = false;
        app.m_simulating = false;
        std::mt19937 rng(std::random_device{}());
        std::uniform_real_distribution<double> pos_dist(-10.0, 10.0);
        app.pin_pos_0 = Vector2F(pos_dist(rng), pos_dist(rng));
        app.pin_pos_1 = Vector2F(pos_dist(rng), pos_dist(rng));
        std::uniform_real_distribution<double> angle_dist(-std::numbers::pi, std::numbers::pi);
        app.pin_angle_0 = angle_dist(rng);
        app.pin_angle_1 = angle_dist(rng);
        app.pins.clear();
        app.setPins();
    }

    def_grad_test_case.input_x = app.x_prev;
    add_tri_energy_test_case.input_x = app.x_prev;
    add_pin_energy_test_case.input_x = app.x_prev;
    add_pin_grad_test_case.input_x = app.x_prev;
    add_pin_hess_test_case.input_x = app.x_prev;
    add_grav_energy_test_case.input_x = app.x_prev;
    add_grav_grad_test_case.input_x = app.x_prev;
    add_grav_hess_test_case.input_x = app.x_prev;
    compute_grad_fd_test_case.input_x = app.x_prev;
    compute_hess_fd_test_case.input_x = app.x_prev;
    add_sys_energy_test_case.input_x = app.x_prev;
    add_sys_grad_test_case.input_x = app.x_prev;
    add_sys_hess_test_case.input_x = app.x_prev;

    add_pin_energy_test_case.pins = app.pins;
    add_pin_grad_test_case.pins = app.pins;
    add_pin_hess_test_case.pins = app.pins;
    compute_grad_fd_test_case.pins = app.pins;
    compute_hess_fd_test_case.pins = app.pins;
    add_sys_energy_test_case.pins = app.pins;
    add_sys_grad_test_case.pins = app.pins;
    add_sys_hess_test_case.pins = app.pins;

    for (size_t i = 0; i < app.elements.size(); ++i) {
        Matrix2F F = app.deformationGradient(app.x_prev, i);
        def_grad_test_case.expected_F.push_back(F);
    }

    for (size_t i = 0; i < app.elements.size(); ++i) {
        double E = 0.0;
        app.addTriangleEnergy(app.x_prev, i, E);
        add_tri_energy_test_case.expected_energies.push_back(E);
    }

    for (size_t i = 0; i < app.pins.size(); ++i) {
        double E = 0.0;
        app.addPinEnergy(app.x_prev, i, E);
        add_pin_energy_test_case.expected_energies.push_back(E);
    }

    VectorXF g = VectorXF::Zero(app.x_prev.size());
    size_t pin_index = 0;
    app.addPinGradient(app.x_prev, pin_index, g);
    add_pin_grad_test_case.expected_pin_index = pin_index;
    add_pin_grad_test_case.expected_gradient = g;

    TripletListF hessian_triplets;
    size_t pin_index_hess = 0;
    app.addPinHessian(app.x_prev, pin_index_hess, hessian_triplets);
    add_pin_hess_test_case.expected_pin_index = pin_index_hess;
    add_pin_hess_test_case.expected_hessians = hessian_triplets;

    double E_grav = 0.0;
    app.addGravityEnergy(app.x_prev, E_grav);
    add_grav_energy_test_case.expected_energies.push_back(E_grav);

    VectorXF g_grav = VectorXF::Zero(app.x_prev.size());
    app.addGravityGradient(app.x_prev, g_grav);
    add_grav_grad_test_case.expected_gradient = g_grav;

    TripletListF hessian_triplets_grav;
    app.addGravityHessian(app.x_prev, hessian_triplets_grav);
    add_grav_hess_test_case.expected_hessians = hessian_triplets_grav;

    VectorXF g_fd = VectorXF::Zero(app.x_prev.size());
    double eps = 1e-8;
    app.computeGradientFD(app.x_prev, g_fd, eps);
    compute_grad_fd_test_case.eps = eps;
    compute_grad_fd_test_case.expected_gradient = g_fd;

    MatrixXF H = MatrixXF::Zero(app.x_prev.size(), app.x_prev.size());
    double eps_hess = 1e-6;
    app.computeHessianFD(app.x_prev, H, eps_hess);
    compute_hess_fd_test_case.eps = eps_hess;
    compute_hess_fd_test_case.expected_H = H;

    double E_sys = 0.0;
    app.addSystemEnergy(app.x_prev, E_sys);
    add_sys_energy_test_case.expected_energies.push_back(E_sys);
    add_sys_energy_test_case.use_gravity = app.use_gravity;

    VectorXF g_sys = VectorXF::Zero(app.x_prev.size());
    app.addSystemGradient(app.x_prev, g_sys);
    add_sys_grad_test_case.expected_gradient = g_sys;

    TripletListF hessian_triplets_sys;
    app.addSystemHessian(app.x_prev, hessian_triplets_sys);
    add_sys_hess_test_case.expected_hessians = hessian_triplets_sys;

    test_cases.push_back(add_tri_energy_test_case);
    test_cases.push_back(def_grad_test_case);
    test_cases.push_back(add_pin_energy_test_case);
    test_cases.push_back(add_pin_grad_test_case);
    test_cases.push_back(add_pin_hess_test_case);
    test_cases.push_back(add_grav_energy_test_case);
    test_cases.push_back(add_grav_grad_test_case);
    test_cases.push_back(add_grav_hess_test_case);
    test_cases.push_back(compute_grad_fd_test_case);
    test_cases.push_back(compute_hess_fd_test_case);
    test_cases.push_back(add_sys_energy_test_case);
    test_cases.push_back(add_sys_grad_test_case);
    test_cases.push_back(add_sys_hess_test_case);

    std::ofstream fout("femTests.dat");
    if (!fout) {
        throw std::runtime_error("Could not open file!");
    }
    for (const auto& test_case : test_cases) {
        test_case.save(fout);
    }

    fout.close();
}


std::string FEMApp::m_check() {

    std::string check_info;

    // generateFEMTests(*this, true);
    
    std::map<std::string, std::vector<TestCase>> test_cases;
    std::ifstream fin(PBS_DATA_DIR "/TestCases.dat");
    if (!fin) {
        return "ERROR: Could not open test file '" PBS_DATA_DIR "/TestCases.dat'!";
    }
    std::string line;
    std::vector<std::string> lines;
    while (std::getline(fin, line)) {
        lines.push_back(line);
        TestCase test_case;
        if (test_case.load(line)) {
            std::string test_name = test_case.test_name;
            test_cases[test_name].push_back(test_case);
        } else {
            check_info += "Failed to load test line: " + line + "\n";
        }
    }
    fin.close();

    for (auto &test : tests) {
        const std::string &test_name = test.first;
        auto it = test_cases.find(test_name);
        if (it != test_cases.end()) {
            FEMApp mockupFEMApp;
            mockupFEMApp.nx = 5;
            mockupFEMApp.ny = 5;
            mockupFEMApp.use_dynamics = false;
            mockupFEMApp.m_simulating = false;
            mockupFEMApp.resetSimulation();
            check_info += test.second(mockupFEMApp, it->second);
        }
    }

    return check_info;
}
