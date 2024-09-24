PYTHON = python3
SCRIPT = convection_diffusion.py
OUTPUT = result.log

all: run

run:
	$(PYTHON) $(SCRIPT) > $(OUTPUT) 2>&1

clean:
	rm -f $(OUTPUT)
