PYTHON = python3
SCRIPT = convection_diffusion.py
OUTPUT = result.log

all: run_outputfile

run_terminal:
	$(PYTHON) $(SCRIPT)

run_outputfile:
	$(PYTHON) $(SCRIPT) > $(OUTPUT) 2>&1

clean:
	rm -f $(OUTPUT)
