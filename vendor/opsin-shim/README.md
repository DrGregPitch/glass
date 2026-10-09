# opsin-shim

`GlassOpsin.java` is a small class compiled into OPSIN's own package so it can run OPSIN's
parse pipeline step by step. For each name it reports the raw token tree (locants, groups,
suffixes, brackets, in name order) and, for every token that created a fragment, the atom
ids of that fragment — the same ids OPSIN's CML uses. `core/names.py` turns that into the
word → atom regions behind name hover, with a vocabulary-based fallback when the sidecar is
unavailable. It also parses with radicals allowed, which is how substituent names (methyl,
phenyl) resolve.

The compiled class is committed (Java 8 bytecode; any JRE 8+ runs it, no JDK needed at
deploy time). To rebuild after editing the source:

```
JAR=$(python -c "import py2opsin,glob,os;print(glob.glob(os.path.join(os.path.dirname(py2opsin.__file__),'opsin*.jar'))[0])")
javac -nowarn -cp "$JAR" -d vendor/opsin-shim vendor/opsin-shim/GlassOpsin.java
```

Protocol: one name per stdin line, one JSON object per stdout line (`core/opsin.py:_ShimWorker`).
