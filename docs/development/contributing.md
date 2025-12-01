# Contributing

Thank you for your interest in contributing to C-GEM Network!

## Getting Started

1. Fork the repository
2. Clone your fork
3. Create a feature branch
4. Make your changes
5. Submit a pull request

## Development Setup

```powershell
# Clone repository
git clone https://github.com/yourusername/CGEM_network.git
cd CGEM_network

# Build
.\scripts\build.bat

# Run tests
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt
```

## Code Guidelines

See [Code Style](code-style.md) for detailed guidelines.

### Summary

- **Language**: ANSI C (C11)
- **Memory**: Allocate at init, free at shutdown
- **Config**: No hardcoded values
- **Documentation**: Comment all functions

## Submitting Changes

### Pull Request Process

1. Update documentation if needed
2. Add tests for new features
3. Ensure build passes without warnings
4. Update CHANGELOG.md
5. Submit PR with clear description

### Commit Messages

```
type(scope): brief description

Detailed explanation if needed.

References: #issue_number
```

Types: `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`

## Reporting Issues

- Use GitHub Issues
- Include: OS, C compiler version, reproduction steps
- Attach minimal test case if possible
