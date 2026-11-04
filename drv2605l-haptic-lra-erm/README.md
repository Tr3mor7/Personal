# drv2605l-haptic-lra-erm: haptic breathing pacer on a DRV2605L driver

An Arduino sketch for a vibration-based breathing pacer. A small motor
vibrates more strongly for 4 s (breathe in), then fades for 6 s (breathe
out), which paces about 6 breaths a minute. The motor is driven through a
**TI DRV2605L haptic driver** over I²C, and one constant chooses between two
kinds of vibration motor.

## Hardware options

### Vibration motor (`USE_ERM` at the top of the sketch)

| | `USE_ERM = false`: **LRA** | `USE_ERM = true`: **ERM** (default) |
|---|---|---|
| Part | MIKROE-6809 (INEED LRA0825BC-0167F) linear resonant actuator | Adafruit #1201 Vibrating Mini Motor Disc (eccentric rotating mass) |
| DRV2605L mode | **closed loop**: auto-resonance tracking, automatic overdrive and braking | **open loop**: no calibration needed |
| Setup at start | auto-calibration (about 1 s); prints its results | none |
| Ratings used | 1.2 V rms at 240 Hz, 1.70 V peak overdrive clamp | 2.0–3.0 V (ESP32) or 2.0–3.3 V (Uno/Nano) |
| Feel | fixed pitch; vibrates out of its face; starts and stops sharply | pitch rises with strength; vibrates in the plane of the disc |
| `f` command | measures the resonant frequency the driver locked onto | sweeps the voltage to find where the motor starts turning |
| Effect library | 6 (LRA) | 2 (ERM library B) |

### Board (chosen automatically when you compile)

| | ESP32 DevKit V1 (ESP32-WROOM-32) | Arduino Uno R3/R4 or Nano |
|---|---|---|
| DRV2605L VIN | **3V3** (the ESP32's pins aren't 5 V tolerant) | 5V |
| SDA / SCL | GPIO21 / GPIO22 | A4 / A5 |
| Start/stop button | on-board **BOOT** button (GPIO0) | push button from pin 2 to GND |
| Buzzer (optional) | GPIO18 → 220 Ω → buzzer + | pin 8 → 220 Ω → buzzer + |
| ERM full-scale voltage | 3.0 V | 3.3 V |

Every board: DRV2605L GND → GND, and DRV2605L IN → GND (the trigger input
isn't used).

Motor leads go to the DRV2605L's **Motor + / −**. The LRA's two leads can go
either way round; the ERM's red lead goes to +.

## Software

* **Arduino IDE**, with the *esp32 by Espressif* core 3.x (board *ESP32 Dev
  Module*) for the ESP32.
* **Adafruit DRV2605 Library** from the Library Manager (it also installs
  Adafruit BusIO).

Open `DRV2605L_Haptic_LRA_ERM/DRV2605L_Haptic_LRA_ERM.ino`, set `USE_ERM`,
upload, and open the Serial Monitor at **115200 baud**. Fix the motor in its
housing before you start: the LRA's calibration fails if the motor is loose.

## Commands (Serial Monitor)

| Key | Action |
|---|---|
| `c` | click (library effect 1, Strong Click) |
| `b` / `h` | buzz at 100% / 50% for 1 s |
| `r` | one slow breath (4 s up, 6 s down) |
| `d` | demo session: chime, 3 breaths, chime |
| `f` | LRA: measure resonance · ERM: start-voltage sweep |
| `x` | re-run calibration (LRA only) |
| `t` | chime on the buzzer |
| `s` | stop |
| `?` | show the menu |

Wait for one command to finish before you send the next. The button (BOOT on
the ESP32) starts the demo session, and pressing it again stops it. On the
ESP32, don't hold BOOT while the board resets, because that puts it into
upload mode.

## How the register values are worked out

`computeRegisters()` calculates the driver settings from the motor's ratings,
using the equations in the DRV2605L datasheet (TI SLOS854), instead of
hard-coding them:

| Register | Equation | Result |
|---|---|---|
| `RATED_VOLTAGE` (0x16), LRA | V<sub>rms</sub> = 20.58 mV × RV / √(1 − (4·t<sub>sample</sub> + 300 µs)·f) | 1.2 V rms at 240 Hz → **0x2F** |
| `OD_CLAMP` (0x17), LRA | V<sub>peak</sub> = 21.22 mV × OD | 1.70 V → **0x50** |
| `OD_CLAMP` (0x17), ERM open loop | V = 21.59 mV × OD | 3.0 V → **0x8B**, 3.3 V → **0x99** |
| `DRIVE_TIME` (0x1B), LRA | half the LRA period = 0.1 ms × DT + 0.5 ms | 2.08 ms → **16** |

For the ERM, strengths of 1–100% are spread between `ERM_V_MIN` (2.0 V) and
full scale, so the motor never stalls at the bottom of a breath. Run `f` to
find your motor's start voltage and set `ERM_V_MIN` just above it.

## AI disclosure

This code was written with the support of Claude Code (Anthropic's AI coding
assistant). The project and its design were directed and reviewed by Gurnoor Hayer.
