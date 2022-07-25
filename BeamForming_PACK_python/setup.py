import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="BeamFormingPyProject", 
    version="0.0.1",
    author="M Yommi",
    author_email="myommi@inti.gob.ar",
    description="Beam forming algorithm using phase-only random masks",
    long_description="Here is presented an iterative algorithm that locally resolves the beam shaping problem through determining high-accuracy random phase-only masks. The algorithm is based on the discrimination of the target field spatial frequencies and a relation between the intensities of the incident and target fields, which is determined from a simple physical model.",
    url="https://www.researchgate.net/publication/337211653_Determining_high-accuracy_random_phase-only_masks_for_complex_local_modulation_of_arbitrary_light_fields",
    packages=['BeamFormingPyProject'],
    python_requires='>=3.8',
)
