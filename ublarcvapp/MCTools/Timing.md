| to &rarr; | electronics time | (_ticks_) | TPC time ticks | trigger time | trigger clock ticks | beam gate time | Optical clock ticks | External clock ticks |
| ------------------------------------------------------------: | --------------------- | -------------------- | ------------------ | ------------------------- | ------------------- | ------------------------- | --------------------- | ---------------------- |
| (unit)                                                        | &micro;s              |                      |                    | &micro;s                  |                     | &micro;s                  |                       |                        |
| @ref DetectorClocksHardwareTrigger         "hardware trigger" | `TriggerTime()`       | `TPCClock()`         |                    |                           | `TriggerClock()`    |                           | `OpticalClock()`      | `ExternalClock()`      |
| @ref DetectorClocksBeamGateOpening         "beam gate point"  | `BeamGateTime()`      |                      |                    |                           |                     |                           |                       |                        |
| @ref DetectorClocksElectronicsTime         "electronics time" |                       |                      | `Time2Tick()`      |                           |                     |                           |                       |                        |
|                                            &nbsp; (_ticks_)   |                       |                      | `TPCTDC2Tick()`    |                           |                     |                           |                       |                        |
| @ref DetectorClocksTPCelectronicsTime      "TPC time"         |                       |                      |                    |                           |                     |                           |                       |                        |
|                                            &nbsp; (_ticks_)   | `TPCTick2Time()`      | `TPCTick2TDC()`      |                    | `TPCTick2TrigTime()`      |                     | `TPCTick2BeamTime()`      |                       |                        |
| @ref DetectorClocksTriggerTime             "trigger time"     |                       |                      |                    |                           |                     |                           |                       |                        |
|                                            &nbsp; (_ticks_)   |                       |                      |                    |                           |                     |                           |                       |                        |
| @ref DetectorClocksOpticalElectronicsTime  "Optical"          |                       |                      |                    |                           |                     |                           |                       |                        |
|                                            &nbsp; (_ticks_)   | `OpticalTick2Time()`  | `OpticalTick2TDC()`  |                    | `OpticalTick2TrigTime()`  |                     | `OpticalTick2BeamTime()`  |                       |                        |
| @ref DetectorClocksExternalElectronicsTime "External"         |                       |                      |                    |                           |                     |                           |                       |                        |
|                                            &nbsp; (_ticks_)   | `ExternalTick2Time()` | `ExternalTick2TDC()` |                    | `ExternalTick2TrigTime()` |                     | `ExternalTick2BeamTime()` |                       |                        |
| @ref DetectorClocksSimulationTime          "simulation time"  | `G4ToElecTime()`      | `TPCG4Time2TDC()`    | `TPCG4Time2Tick()` |                           |                     |                           | `OpticalG4Time2TDC()` | `ExternalG4Time2TDC()` |