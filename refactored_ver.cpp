#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <memory>

class Vector3D
{
public:
    double x, y, z;

    Vector3D() : x(0.0), y(0.0), z(0.0) {}

    Vector3D(double xVal, double yVal, double zVal) : x(xVal), y(yVal), z(zVal) {}

    Vector3D operator+(const Vector3D &other) const
    {
        return Vector3D(x + other.x, y + other.y, z + other.z);
    }

    Vector3D operator-(const Vector3D &other) const
    {
        return Vector3D(x - other.x, y - other.y, z - other.z);
    }

    double dotProduct(const Vector3D &other) const
    {
        return x * other.x + y * other.y + z * other.z;
    }

    Vector3D crossProduct(const Vector3D &other) const
    {
        return Vector3D(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
    }

    Vector3D unitVector() const
    {
        double magnitude = std::sqrt(x * x + y * y + z * z);
        if (magnitude != 0.0)
        {
            return Vector3D(x / magnitude, y / magnitude, z / magnitude);
        }
        else
        {
            // Handle the case where the vector is zero
            return *this; // Return the zero vector itself
        }
    }

    Vector3D scalarMultiply(double scalar) const
    {
        return Vector3D(x * scalar, y * scalar, z * scalar);
    }
};
class Matrix3x3
{
public:
    std::vector<std::vector<double>> data;

    Matrix3x3() : data(3, std::vector<double>(3, 0.0)) {}

    Matrix3x3(const std::vector<std::vector<double>> &inputData) : data(inputData) {}

    Vector3D operator*(const Vector3D &vec) const
    {
        if (data.size() != 3 || data[0].size() != 3)
        {
            throw std::invalid_argument("Matrix dimensions are invalid.");
        }

        double resultX = data[0][0] * vec.x + data[0][1] * vec.y + data[0][2] * vec.z;
        double resultY = data[1][0] * vec.x + data[1][1] * vec.y + data[1][2] * vec.z;
        double resultZ = data[2][0] * vec.x + data[2][1] * vec.y + data[2][2] * vec.z;

        return Vector3D(resultX, resultY, resultZ);
    }

    Vector3D unitVector() const
    {
        double rowMagnitudes[3] = {
            std::sqrt(data[0][0] * data[0][0] + data[0][1] * data[0][1] + data[0][2] * data[0][2]),
            std::sqrt(data[1][0] * data[1][0] + data[1][1] * data[1][1] + data[1][2] * data[1][2]),
            std::sqrt(data[2][0] * data[2][0] + data[2][1] * data[2][1] + data[2][2] * data[2][2])};

        if (rowMagnitudes[0] != 0.0 && rowMagnitudes[1] != 0.0 && rowMagnitudes[2] != 0.0)
        {
            return Vector3D(data[0][0] / rowMagnitudes[0], data[1][0] / rowMagnitudes[1], data[2][0] / rowMagnitudes[2]);
        }
        else
        {
            throw std::invalid_argument("Cannot normalize zero vector(s).");
        }
    }
};

class PhysicsCalculations
{
public:
    static constexpr double GRAVITY = 9.81;            // Acceleration due to gravity in m/s^2
    static constexpr double SPECIFIC_IMPULSE = 2500.0; // Specific impulse in seconds
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double LUNAR_GRAVITY = 1.625; // Gravity on the Moon's surface in m/s^2
    static constexpr double TTFSCALE = 1.0;        // Define TTFSCALE value
    static constexpr double TREDESCL = -0.08;      // Define TREDESCL value
    static constexpr double NEGMAX = -1.0;         // Define NEGMAX value
    static constexpr double POSMAX = 1.0;          // Define POSMAX value

    static double calculateVelocity(double initialVelocity, double acceleration, double time)
    {
        return initialVelocity + acceleration * time;
    }

    static double calculateAltitudeRate(double velocityY)
    {
        return velocityY;
    }

    static double calculateHorizontalDistance(double velocity, double time)
    {
        return velocity * time;
    }

    static double calculateFuelConsumptionRate(double thrust)
    {
        if (SPECIFIC_IMPULSE != 0.0)
        {
            return thrust / SPECIFIC_IMPULSE;
        }
        else
        {
            throw std::invalid_argument("Specific impulse cannot be zero.");
        }
    }

    static double calculateGravityForce(double mass, bool isMoon = false)
    {
        double gravity = isMoon ? LUNAR_GRAVITY : GRAVITY;
        return gravity * mass;
    }

    static double calculateDragForce(double dragCoefficient, double airDensity, double velocity)
    {
        return 0.5 * dragCoefficient * airDensity * velocity * velocity;
    }

    static double calculateTrajectoryAngle(double horizontalVelocity, double verticalVelocity)
    {
        if (horizontalVelocity != 0.0)
        {
            return std::atan2(verticalVelocity, horizontalVelocity);
        }
        else
        {
            throw std::invalid_argument("Horizontal velocity cannot be zero for angle calculation.");
        }
    }

    static double calculateRange(double horizontalVelocity, double time)
    {
        return horizontalVelocity * time;
    }
};

class GuidanceSystem
{
private:
    int wchPhase;
    double ttfOver8;
    double gain;
    Vector3D land;
    Vector3D r;
    Vector3D cg;
    Vector3D vgu;
    Vector3D angTerm;
    Matrix3x3 transformationMatrix;

public:
    GuidanceSystem() : wchPhase(0), ttfOver8(0.0), gain(1.0)
    {
        land = Vector3D();
        r = Vector3D();
        cg = Vector3D();
        vgu = Vector3D();
        angTerm = Vector3D();
        transformationMatrix = Matrix3x3();
    }

    void setWchPhase(int phase)
    {
        wchPhase = phase;
    }

    void setTtfOver8(double value)
    {
        ttfOver8 = value;
    }

    void setGain(double value)
    {
        gain = value;
    }

    void setLand(const Vector3D &vec)
    {
        land = vec;
    }

    void setR(const Vector3D &vec)
    {
        r = vec;
    }

    Vector3D unitVector() const
    {
        // Calculate the magnitude of each row vector
        double rowMagnitudes[3] = {
            std::sqrt(r.x * r.x + r.y * r.y + r.z * r.z),
            std::sqrt(cg.x * cg.x + cg.y * cg.y + cg.z * cg.z),
            std::sqrt(vgu.x * vgu.x + vgu.y * vgu.y + vgu.z * vgu.z)};

        // Check if any magnitude is zero
        if (rowMagnitudes[0] != 0.0 && rowMagnitudes[1] != 0.0 && rowMagnitudes[2] != 0.0)
        {
            // Return the normalized row vectors as columns of the unit matrix
            return Vector3D(r.x / rowMagnitudes[0], cg.x / rowMagnitudes[1], vgu.x / rowMagnitudes[2]);
        }
        else
        {
            // Handle the case where any row vector is zero
            return Vector3D(); // Return the zero vector
        }
    }

    void calculateAngTerm()
    {
        // Calculate angTerm using vector operations
        Vector3D vguUnit = vgu.unitVector();                                 // Get the unit vector of vgu
        Vector3D crossProd = r.crossProduct(transformationMatrix * vguUnit); // Compute cross product
        angTerm = crossProd.scalarMultiply(PhysicsCalculations::PI);         // Multiply by PI to get angTerm
    }

    void computeCg()
    {
        cg = transformationMatrix * (r - land);
    }

    void computeTtfOver8()
    {
        ttfOver8 = ttfOver8 / 8.0;
    }

    void computeMainGuidanceEquation()
    {
        // Main guidance equation calculations
    }

    void erectTransformationMatrix()
    {
        // Erect guidance-stable member transformation matrix
        // Compute and update transformationMatrix
    }

    void prepareForExit()
    {
        // Prepare for system exit based on phase and conditions
    }

    void exitLandingGuidance()
    {
        // Handle system exit during landing guidance phases
    }

    void displayInformation()
    {
        // Display system information based on phase and parameters
    }

    void handleAlarm()
    {
        // Alarm handling routine for specific conditions
    }

    void handleOverflow()
    {
        // Overflow handling routine
    }

    void controlThrottle()
    {
        // Throttle control logic
    }

    void findCduw()
    {
        // Find CDUW (Control Display Unit Word) based on specific calculations
    }

    void integrateLunarLandingEquations()
    {
        // Lunar landing constants
        const double TSCALINV = 1.0 / PhysicsCalculations::TTFSCALE;
        const double DEC103 = -103.0;
        const double DEC99 = 99.0;
        const double TREDESCL = PhysicsCalculations::TREDESCL;
        const double NEGMAX = PhysicsCalculations::NEGMAX; // Define NEGMAX value
        const double POSMAX = PhysicsCalculations::POSMAX; // Define POSMAX value

        // Example calculation using provided assembly code comments
        double result = ttfOver8 / 8.0 * TREDESCL - DEC103 + NEGMAX + POSMAX;
        std::cout << "Result: " << result << std::endl;
    }
};

int main()
{
    GuidanceSystem guidanceSystem;

    guidanceSystem.setWchPhase(1);
    guidanceSystem.setTtfOver8(2.5);
    guidanceSystem.setGain(0.5);

    guidanceSystem.calculateAngTerm();
    guidanceSystem.computeCg();
    guidanceSystem.computeTtfOver8();
    guidanceSystem.computeMainGuidanceEquation();
    guidanceSystem.erectTransformationMatrix();
    guidanceSystem.prepareForExit();
    guidanceSystem.exitLandingGuidance();
    guidanceSystem.displayInformation();
    guidanceSystem.handleAlarm();
    guidanceSystem.handleOverflow();
    guidanceSystem.controlThrottle();
    guidanceSystem.findCduw();

    guidanceSystem.integrateLunarLandingEquations();

    double initialVelocity = 10.0; // m/s
    double acceleration = -2.0;    // m/s^2
    double time = 5.0;             // seconds
    double velocity = PhysicsCalculations::calculateVelocity(initialVelocity, acceleration, time);
    std::cout << "Velocity after " << time << " seconds: " << velocity << " m/s" << std::endl;

    return 0;
}
