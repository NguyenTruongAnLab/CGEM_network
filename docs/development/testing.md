# Testing

## Test Cases

### Available Cases

| Case | Purpose | Complexity |
|------|---------|------------|
| `Tien_River` | Basic multi-branch | Simple |
| `Mekong_Delta_Full` | Full validation | Complex |
| `SaigonDongNai` | Tributary handling | Medium |

### Running Tests

```powershell
# Run all test cases
.\scripts\build-and-run.ps1 -r Tien_River
.\scripts\build-and-run.ps1 -r Mekong_Delta_Full
.\scripts\build-and-run.ps1 -r SaigonDongNai
```

## Validation Checks

### 1. Mass Conservation

```python
# Check mass balance at junctions
def check_mass_balance(branches, junction_node):
    Q_in = sum(b['Q'] for b in branches if b['direction'] == 'in')
    Q_out = sum(b['Q'] for b in branches if b['direction'] == 'out')
    assert abs(Q_in - Q_out) < 0.01 * Q_in, "Mass balance violation"
```

### 2. Numerical Stability

- No NaN values in output
- No negative concentrations
- Bounded velocities

### 3. Physical Consistency

- Salinity between 0-35 PSU
- O₂ positive
- pH between 6-9

## Adding New Tests

1. Create case directory under `INPUT/Cases/`
2. Add configuration files
3. Document expected results
4. Add to CI workflow
