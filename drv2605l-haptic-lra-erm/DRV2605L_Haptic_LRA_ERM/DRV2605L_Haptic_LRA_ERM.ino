/*
  DRV2605L_Haptic_LRA_ERM.ino
  Haptic breathing-pacer prototype: TI DRV2605L haptic driver with a choice of
  vibration motor (LRA in closed loop, or ERM in open loop) on ESP32 or Arduino.
  (Originally "Version A" of the BMED 2310 Team AlphaX prototype: the version
  that uses a DRV2605L driver chip.)

  Drives one of two motors through a TI DRV2605L haptic driver (Adafruit
  DRV2605L breakout, #2305) over I2C. Pick the motor with USE_ERM below.

  USE_ERM = false: MIKROE-6809 linear resonant actuator (INEED LRA0825BC-0167F)
    The driver finds and tracks the LRA's resonance by itself and adds automatic
    overdrive and braking, so the motor starts and stops quickly.
    Ratings (LRA0825BC-0167F datasheet): 1.2 +/- 0.05 Vrms sine at 240 +/- 10 Hz,
    0.1-1.25 Vrms operating range, <= 90 mA rms, 20 +/- 3 ohm coil, vibrates
    along the Z axis (out of the face).

  USE_ERM = true: Adafruit Vibrating Mini Motor Disc (#1201), an eccentric
    rotating mass (ERM) motor: a small off-center weight spun by a DC motor.
    Ratings (Adafruit product page): 10 mm x 2.7 mm, 0.9 g, rated 2.5-3.8 V,
    vibrates from about 2 V to 5 V, 60 mA at 3 V, 11,000 RPM at 5 V.
    The driver runs it open loop (no calibration needed). Unlike the LRA, its
    buzz gets higher in pitch as well as stronger as the voltage rises, and it
    shakes in the plane of the disc rather than out of its face.

  Board: ESP32 DevKit V1 (ESP32-WROOM-32), programmed in the Arduino IDE with
  the "esp32 by Espressif" core 3.x, board "ESP32 Dev Module". The sketch still
  runs on an Uno R3/R4 or Nano; the pins switch automatically (see Pins below).

  Wiring (ESP32 DevKit V1):
    DRV2605L VIN   -> 3V3  (not VIN/5V: the breakout's I2C pull-ups go to its
                            VIN, and the ESP32's pins are not 5 V tolerant)
    DRV2605L GND   -> GND
    DRV2605L SDA   -> GPIO21 (D21)
    DRV2605L SCL   -> GPIO22 (D22)
    DRV2605L IN    -> GND  (trigger input, not used here)
    Motor + / -    -> LRA: the two leads, either way round (it is AC-driven)
                      ERM: red lead to Motor +, blue lead to Motor -
    GPIO18 (D18)   -> 220 ohm resistor -> buzzer +   (optional passive buzzer)
    Buzzer -       -> GND
    BOOT button (GPIO0, on the board) starts and stops the demo; no wiring needed.
  On an Uno or Nano instead: VIN -> 5V, SDA -> A4, SCL -> A5, buzzer on pin 8,
  and a push button from pin 2 to GND.

  Library: "Adafruit DRV2605 Library" from the Arduino Library Manager
  (it also installs Adafruit BusIO).

  Use: fix the motor in its housing first, then upload and open the Serial
  Monitor at 115200 baud. Send one letter:
    c = click            b = buzz 100%, 1 s      h = buzz 50%, 1 s
    r = one slow breath (4 s up, 6 s down)       d = demo session
    f = LRA: measure resonance / ERM: start-voltage sweep
    x = recalibrate (LRA only)
    t = chime            s = stop                ? = this list
  Wait for one command to finish before sending the next.
  The BOOT button (pin 2 button on an Uno) starts the demo session; press it
  again to stop. Press it only after the board has started: holding BOOT while
  the board resets puts it in upload mode.
*/

#include <Wire.h>
#include <Adafruit_DRV2605.h>

Adafruit_DRV2605 drv;

// ---------------- Motor choice ----------------
// true  = Adafruit Vibrating Mini Motor Disc (#1201), ERM
// false = MIKROE-6809 LRA
const bool USE_ERM = true;

// ---------------- Pins ----------------
// SDA and SCL are the board's defaults (GPIO21 / GPIO22 on the ESP32 DevKit V1).
#if defined(ARDUINO_ARCH_ESP32)
const int BUZZER_PIN = 18;  // optional passive buzzer, through 220 ohm
const int BUTTON_PIN = 0;   // the board's BOOT button (GPIO0, already pulled up)
#else
const int BUZZER_PIN = 8;   // Uno / Nano: optional passive buzzer, through 220 ohm
const int BUTTON_PIN = 2;   // Uno / Nano: optional push button to GND
#endif

// ---------------- Settings ----------------
const int DEMO_BREATHS = 3;  // breaths in the demo session, 10 s each

// ---------------- LRA values (LRA0825BC-0167F datasheet) ----------------
const float LRA_VRMS  = 1.20;   // rated drive voltage, V rms (max 1.25 V rms)
const float LRA_FREQ  = 240.0;  // resonant frequency, Hz
const float LRA_VPEAK = 1.70;   // overdrive clamp, V peak. No overdrive rating is
                                // given, so the peak of the rated sine is used
                                // (1.2 V rms x 1.414), per TI app note SLOA189.

// DRV2605L back-EMF sample time used in the LRA rated-voltage equation
const float SAMPLE_TIME_S = 300e-6;  // SAMPLE_TIME = 3 (300 us, the default)

// ---------------- ERM values (Adafruit #1201 product page) ----------------
#if defined(ARDUINO_ARCH_ESP32)
const float ERM_V_FULL = 3.0;  // V at 100% strength. The DRV2605L runs from the ESP32's
                               // 3.3 V, so its output cannot go much past 3.2 V; 3.0 V
                               // is also the chip's default full scale.
#else
const float ERM_V_FULL = 3.3;  // V at 100% strength with the DRV2605L on 5 V.
#endif                         // Both are inside the motor's 2.5-3.8 V rating.
const float ERM_V_MIN  = 2.0;  // V at the lowest non-zero strength. The motor stalls
                               // below roughly 2 V, so strengths of 1-100% are spread
                               // over ERM_V_MIN to ERM_V_FULL (50% = halfway between).
                               // Run 'f' to find your motor's start voltage and set
                               // this just above it.
const uint8_t ERM_LIBRARY = 2; // DRV2605L ERM library B: 3 V rated, 3 V overdrive,
                               // 40-60 ms rise (closest to this motor). Used by 'c'.

uint8_t ratedVoltageReg;  // register 0x16 (LRA only; ignored in ERM open loop)
uint8_t odClampReg;       // register 0x17
uint8_t driveTimeReg;     // register 0x1B, bits 4:0 (LRA only)

// Round a calculated value into one register byte.
uint8_t toByte(float x) {
  if (x < 0) x = 0;
  if (x > 255) x = 255;
  return (uint8_t)(x + 0.5);
}

// Register values from the DRV2605L datasheet (SLOS854).
void computeRegisters() {
  if (USE_ERM) {
    // Eq. 6 (ERM open loop): full-scale V = 21.59 mV x OD_CLAMP
    odClampReg = toByte(ERM_V_FULL / 0.02159);         // 3.0 V -> 139 (0x8B); 3.3 V -> 153 (0x99)
    ratedVoltageReg = 0;                               // not used in open loop
    driveTimeReg = 0;                                  // not used in open loop
    return;
  }

  // Eq. 5: Vrms = 20.58 mV x RATED_VOLTAGE / sqrt(1 - (4 x t_sample + 300 us) x f_LRA)
  float k = sqrt(1.0 - (4.0 * SAMPLE_TIME_S + 300e-6) * LRA_FREQ);
  ratedVoltageReg = toByte(LRA_VRMS * k / 0.02058);  // 1.2 V at 240 Hz -> 47 (0x2F)

  // Eq. 9: Vpeak = 21.22 mV x OD_CLAMP
  odClampReg = toByte(LRA_VPEAK / 0.02122);          // 1.70 V -> 80 (0x50)

  // Drive time should start at about half the LRA period:
  // drive time (ms) = DRIVE_TIME x 0.1 ms + 0.5 ms
  float halfPeriodMs = 500.0 / LRA_FREQ;             // 2.08 ms at 240 Hz
  driveTimeReg = toByte((halfPeriodMs - 0.5) / 0.1); // -> 16 (2.1 ms)
  if (driveTimeReg > 31) driveTimeReg = 31;
}

// ---------------- Button and stop handling ----------------
bool buttonDown = false;     // debounced button state
bool lastReading = false;
unsigned long lastChangeMs = 0;

// True once each time the button is pressed (20 ms debounce).
bool buttonPressed() {
  bool reading = (digitalRead(BUTTON_PIN) == LOW);
  if (reading != lastReading) {
    lastReading = reading;
    lastChangeMs = millis();
  }
  if (millis() - lastChangeMs > 20 && reading != buttonDown) {
    buttonDown = reading;
    return buttonDown;       // report presses, not releases
  }
  return false;
}

// Checked while a pattern runs: 's' in the Serial Monitor or a button press stops it.
bool stopRequested() {
  while (Serial.available()) {
    if (Serial.read() == 's') return true;
  }
  return buttonPressed();
}

// ---------------- Driver control ----------------
// Low-power standby (about 4 uA). Registers and calibration are kept.
void standby() {
  drv.writeRegister8(DRV2605_REG_MODE, 0x40);
}

// Stop any vibration, let the motor brake (LRA) or spin down (ERM), then standby.
void stopMotor() {
  drv.setRealtimeValue(0);
  drv.stop();
  delay(60);
  standby();
}

// ERM setup: open loop, so no calibration is needed.
void configureErm() {
  drv.writeRegister8(DRV2605_REG_MODE, 0x00);       // leave standby
  // Feedback: ERM mode (bit 7 = 0); the other bits are the chip defaults and
  // only matter in closed loop.
  drv.writeRegister8(DRV2605_REG_FEEDBACK, 0x36);
  drv.writeRegister8(DRV2605_REG_CLAMPV, odClampReg);
  // Control3: ERM_OPEN_LOOP = 1, signed RTP data (0x7F = 100%)
  drv.writeRegister8(DRV2605_REG_CONTROL3, 0xA0);
  drv.selectLibrary(ERM_LIBRARY);
  standby();
}

// LRA auto-calibration (DRV2605L datasheet section 8.5.6). Returns true if it passed.
bool calibrate() {
  drv.writeRegister8(DRV2605_REG_MODE, DRV2605_MODE_AUTOCAL);  // leave standby, calibration mode
  drv.writeRegister8(DRV2605_REG_RATEDV, ratedVoltageReg);
  drv.writeRegister8(DRV2605_REG_CLAMPV, odClampReg);
  // Feedback: LRA mode | brake factor 3x | loop gain high | back-EMF gain 15x (set by calibration)
  drv.writeRegister8(DRV2605_REG_FEEDBACK, 0x80 | (2 << 4) | (2 << 2) | 2);  // 0xAA
  // Control1: startup boost on | DRIVE_TIME
  drv.writeRegister8(DRV2605_REG_CONTROL1, 0x80 | driveTimeReg);             // 0x90
  // Control2: bidirectional input, brake stabilizer, 300 us sample time, blanking/IDISS = 1
  drv.writeRegister8(DRV2605_REG_CONTROL2, 0xF5);
  // Control3: closed loop with auto-resonance (LRA_OPEN_LOOP = 0), signed RTP data
  drv.writeRegister8(DRV2605_REG_CONTROL3, 0xA0);
  // Control4: calibration time 1000-1200 ms. Bit 0 is OTP_PROGRAM and must stay 0.
  drv.writeRegister8(DRV2605_REG_CONTROL4, 0x30);

  drv.go();
  unsigned long t0 = millis();
  while (drv.readRegister8(DRV2605_REG_GO) & 0x01) {  // GO clears when calibration ends
    if (millis() - t0 > 3000) return false;
    delay(20);
  }
  uint8_t status = drv.readRegister8(DRV2605_REG_STATUS);
  return (status & 0x08) == 0;  // DIAG_RESULT (bit 3): 0 = passed
}

void runCalibration() {
  if (USE_ERM) {
    Serial.println(F("Not needed: the ERM runs open loop."));
    return;
  }
  Serial.println(F("Calibrating (about 1 s)..."));
  if (calibrate()) {
    Serial.println(F("Calibration passed."));
    Serial.print(F("  A_CAL_COMP = 0x"));
    Serial.println(drv.readRegister8(DRV2605_REG_AUTOCALCOMP), HEX);
    Serial.print(F("  A_CAL_BEMF = 0x"));
    Serial.println(drv.readRegister8(DRV2605_REG_AUTOCALEMP), HEX);
    Serial.print(F("  BEMF_GAIN  = "));
    Serial.println(drv.readRegister8(DRV2605_REG_FEEDBACK) & 0x03);
  } else {
    Serial.println(F("Calibration FAILED. Check that the LRA is fixed in its housing, the"));
    Serial.println(F("leads are soldered to Motor +/-, and the coil reads about 20 ohm."));
  }
  drv.selectLibrary(6);  // LRA effect library (closed loop)
  standby();
}

// Drive voltage -> RTP value (signed format: 0x7F = full scale)
uint8_t voltsToRtp(float v, float fullScale) {
  if (v <= 0) return 0;
  if (v > fullScale) v = fullScale;
  return (uint8_t)(v / fullScale * 127.0 + 0.5);
}

// Strength in % -> RTP value.
//   LRA: % of the rated voltage (0x7F = 100%).
//   ERM: 1-100% is spread over ERM_V_MIN to ERM_V_FULL so the motor never stalls;
//        0% is off.
uint8_t pctToRtp(float pct) {
  if (pct <= 0) return 0;
  if (pct > 100) pct = 100;
  if (USE_ERM) {
    return voltsToRtp(ERM_V_MIN + (ERM_V_FULL - ERM_V_MIN) * pct / 100.0, ERM_V_FULL);
  }
  return (uint8_t)(pct * 127.0 / 100.0 + 0.5);
}

// Play one effect from the selected library, wait for it to finish, then standby.
void playEffect(uint8_t effect) {
  drv.setMode(DRV2605_MODE_INTTRIG);  // also clears standby
  drv.setWaveform(0, effect);
  drv.setWaveform(1, 0);              // end of sequence
  drv.go();
  unsigned long t0 = millis();
  while ((drv.readRegister8(DRV2605_REG_GO) & 0x01) && millis() - t0 < 2000) {
    delay(5);
  }
  standby();
}

// Change the strength in a straight line from startPct to endPct over ms
// milliseconds (real-time playback). Leaves the motor running.
// Returns false if a stop was requested.
bool ramp(unsigned long ms, float startPct, float endPct) {
  drv.setRealtimeValue(pctToRtp(startPct));
  drv.setMode(DRV2605_MODE_REALTIME);
  unsigned long startMs = millis();
  while (millis() - startMs < ms) {
    float pct = startPct + (endPct - startPct) * (millis() - startMs) / (float)ms;
    drv.setRealtimeValue(pctToRtp(pct));
    if (stopRequested()) return false;
    delay(20);
  }
  return true;
}

// Steady vibration at one strength.
bool buzz(float pct, unsigned long ms) {
  bool done = ramp(ms, pct, pct);
  stopMotor();
  return done;
}

// One slow breath: the vibration rises for 4 s (breathe in) and fades for 6 s
// (breathe out). Ten seconds a breath paces about 6 breaths a minute.
bool breath() {
  bool done = ramp(4000, 10, 100) && ramp(6000, 100, 10);
  stopMotor();
  return done;
}

// LRA only: read the resonance the driver has locked onto (register 0x22).
// The value is only valid while the motor is being driven.
void measureResonance() {
  drv.setRealtimeValue(pctToRtp(100));
  drv.setMode(DRV2605_MODE_REALTIME);
  delay(300);  // give the auto-resonance engine time to lock on
  uint8_t period = drv.readRegister8(DRV2605_REG_LRARESON);
  stopMotor();

  if (period == 0) {
    Serial.println(F("No reading - check the motor connection and mounting."));
    return;
  }
  float periodUs = period * 98.46;  // datasheet: period (us) = LRA_PERIOD x 98.46 us
  Serial.print(F("Resonance: "));
  Serial.print(1e6 / periodUs, 1);
  Serial.println(F(" Hz (datasheet: 230-250 Hz)."));
}

// ERM only: raise the voltage from rest in 0.1 V steps (0.7 s each) and print
// each step. The first step where the motor vibrates is its start voltage;
// set ERM_V_MIN just above it so the bottom of a breath never stalls.
void findStartVoltage() {
  Serial.println(F("Start-voltage sweep: note the first step that vibrates (s stops)"));
  drv.setRealtimeValue(0);
  drv.setMode(DRV2605_MODE_REALTIME);
  for (int mv = 1000; mv <= (int)(ERM_V_FULL * 1000 + 0.5); mv += 100) {
    float v = mv / 1000.0;
    drv.setRealtimeValue(voltsToRtp(v, ERM_V_FULL));
    Serial.print(v, 1);
    Serial.println(F(" V"));
    unsigned long t0 = millis();
    while (millis() - t0 < 700) {
      if (stopRequested()) {
        stopMotor();
        Serial.println(F("Stopped"));
        return;
      }
      delay(20);
    }
  }
  stopMotor();
  Serial.println(F("Sweep done. Set ERM_V_MIN just above the first step that vibrated."));
}

// ---------------- Sound ----------------
// Two-note chime on the buzzer: rising = start, falling = end.
void chime(bool rising) {
  const int lowNote = 1568, highNote = 2093;  // G6 and C7
  tone(BUZZER_PIN, rising ? lowNote : highNote, 120);
  delay(160);
  tone(BUZZER_PIN, rising ? highNote : lowNote, 180);
  delay(220);
}

// ---------------- Commands ----------------
// Demo session: chime, DEMO_BREATHS slow breaths, chime.
void demoSession() {
  Serial.print(F("Demo session: "));
  Serial.print(DEMO_BREATHS);
  Serial.println(F(" slow breaths (s or the button stops it)"));
  chime(true);
  for (int i = 0; i < DEMO_BREATHS; i++) {
    if (!breath()) {
      chime(false);
      Serial.println(F("Stopped"));
      return;
    }
  }
  chime(false);
  Serial.println(F("Demo done"));
}

void printMenu() {
  Serial.println(F("Send: c = click, b = buzz 100%, h = buzz 50%, r = breath, d = demo,"));
  if (USE_ERM) {
    Serial.println(F("      f = start-voltage sweep, t = chime, s = stop, ? = menu"));
  } else {
    Serial.println(F("      f = resonance, x = recalibrate, t = chime, s = stop, ? = menu"));
  }
}

void setup() {
  pinMode(BUTTON_PIN, INPUT_PULLUP);
  Serial.begin(115200);
  while (!Serial && millis() < 2000) {
  }  // wait briefly for boards with native USB

  if (USE_ERM) {
    Serial.println(F("DRV2605L haptic driver - ERM (Adafruit #1201), open loop"));
  } else {
    Serial.println(F("DRV2605L haptic driver - LRA (MIKROE-6809), closed loop"));
  }
  if (!drv.begin()) {
    Serial.println(F("DRV2605L not found. Check VIN, GND, SDA and SCL."));
    while (true) {
    }
  }

  computeRegisters();
  if (USE_ERM) {
    drv.useERM();
    configureErm();
    Serial.print(F("OD_CLAMP = 0x"));
    Serial.print(odClampReg, HEX);
    Serial.print(F(" ("));
    Serial.print(ERM_V_FULL, 1);
    Serial.print(F(" V full scale), library "));
    Serial.println(ERM_LIBRARY);
  } else {
    drv.useLRA();
    Serial.print(F("RATED_VOLTAGE = 0x"));
    Serial.print(ratedVoltageReg, HEX);
    Serial.print(F(", OD_CLAMP = 0x"));
    Serial.print(odClampReg, HEX);
    Serial.print(F(", DRIVE_TIME = "));
    Serial.println(driveTimeReg);
    runCalibration();
  }
  printMenu();
}

void loop() {
  if (buttonPressed()) demoSession();

  if (!Serial.available()) return;
  char cmd = Serial.read();

  switch (cmd) {
    case 'c':
      Serial.println(F("Click"));
      playEffect(1);  // effect 1 = Strong Click, 100%
      break;
    case 'b':
      Serial.println(F("Buzz 100%, 1 s"));
      buzz(100, 1000);
      break;
    case 'h':
      Serial.println(F("Buzz 50%, 1 s"));
      buzz(50, 1000);
      break;
    case 'r':
      Serial.println(F("One breath: 4 s up, 6 s down"));
      breath();
      break;
    case 'd':
      demoSession();
      break;
    case 'f':
      if (USE_ERM) {
        findStartVoltage();
      } else {
        measureResonance();
      }
      break;
    case 'x':
      runCalibration();
      break;
    case 't':
      Serial.println(F("Chime"));
      chime(true);
      break;
    case 's':
      stopMotor();
      Serial.println(F("Stopped"));
      break;
    case '?':
      printMenu();
      break;
    default:
      break;  // ignore line endings and other characters
  }
}
