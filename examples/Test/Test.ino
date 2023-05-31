#include <Quadrature.h>

void setup() {
    Serial.begin(9600);
    
    Serial.println("Calculating integrals on [0, PI]");
    for (int coeff = 0; coeff < 8; coeff++) {
        float result = integrate([coeff](double x) { return sin(2 * coeff + 1) * x; }, 0, PI, 1e-4);
        Serial.print("sin(");
        Serial.print(2 * coeff + 1);
        Serial.print("x): ");
        Serial.println(result, 3);
    }

    Serial.println("Calculating integral of exp(x) on [0, 1]");
    for (int exponent = -3; exponent >= -8; exponent--) {
        float result = integrate(exp, 0, 1, pow(10, exponent));
        Serial.print("precision = 1e");
        Serial.print(exponent);
        Serial.print(": ");
        Serial.println(result, 8);
    }

    Serial.println("Finished!");
}

void loop() {

}