# GAPIT AI User Manual

## What is GAPIT AI

GAPIT AI is an expert assistant skill that runs inside the TRAE AI coding environment. It wraps GAPIT, the Genome Association and Prediction Integrated Tool, into a conversational helper that can design, explain, and generate complete genome-wide association study (GWAS) and genomic selection (GS) analyses. You describe your study in plain English, and the skill returns a working R pipeline instead of forcing you to write every line of code yourself.

The product is also called GAPIT5. The skill folder inside the package is named `gapit-expert`, so you may see all three names used interchangeably. They refer to the same thing.

The skill covers the whole workflow from start to finish: quality control of genotype data, phenotype simulation, population structure, GWAS under multiple models (GLM, MLM, CMLM, MLMM, SUPER, FarmCPU, and BLINK), genomic selection under several methods (gBLUP, cBLUP, sBLUP, GAGBLUP, MAS, BayesB, and rrBLUP), statistical power and false discovery rate, and export of results for downstream use.

## What you need before you start

Three things are required. You need a computer running TRAE, because the skill is loaded and invoked through the TRAE interface. You need a working R installation, because GAPIT runs as R code that you execute locally. And you need the GAPIT AI package itself, which you obtain from the maintainers and place into a TRAE workspace.

Make sure R is available from a terminal, meaning the command `Rscript` works. A few helper packages should also be installed: `BGLR` and `rrBLUP` support the BayesB and rrBLUP genomic selection methods, and `data.table` speeds up file handling. Everything else that makes up GAPIT is loaded from a bundled function file that ships inside the package, so you do not need to install GAPIT itself as a separate R package.

## Downloading and installing TRAE

TRAE is an AI-native development environment. It ships in a few forms. TraeCode is the desktop IDE you will use to open a project and run the GAPIT AI skill. TraeWork is a web, desktop, and mobile workspace. For this workflow you want TraeCode.

Download TraeCode from the official website at `https://www.trae.ai/` and click the Download TRAE button. The official documentation lives at `https://docs.trae.ai/`. If you are based in mainland China, the equivalent site is `https://www.trae.cn/`. After downloading, run the installer and sign in with your account.

You do not need a paid plan to load and try the skill, but your plan's model allowance determines how much analysis you can run. That is explained in the token cost section later in this manual.

## Getting the GAPIT AI package

The skill is distributed as a zip archive named `gapit-ai-v1.0-beta.1.1.zip`. Download it and unzip it into a folder on your computer. The archive contains the skill itself under `.trae/skills/gapit-expert/`, together with the data resources it needs: the GAPIT function file, the two mode-code templates, demo data, and the public question bank.

The single most important file inside is `data/gapit-resources/code/gapit_functions.txt`. This one file holds the entire GAPIT function library. The generated R scripts load it directly with `source()` rather than calling `library(GAPIT)`. That design removes the most common installation problem users reported: R being unable to find the GAPIT package.

## Loading GAPIT AI in TRAE

The skill is a workspace skill. A workspace skill is simply a folder named after the skill, placed inside a `.trae/skills/` directory, which TRAE discovers automatically when you open the containing project.

Open TraeCode and open the folder that holds the unzipped package as your project workspace. If the file `.trae/skills/gapit-expert/SKILL.md` exists in that workspace, TRAE registers the skill automatically. Then mention the skill by name in your conversation and ask your question. For example: "Use the GAPIT AI skill to check my BLINK parameter settings."

If you want the skill available in every project rather than one workspace, copy the skill folder to your user-level skills directory. On macOS and Linux this is `~/.trae/skills/`. The package includes an installer script named `install_gapit_expert.sh` that performs this copy for you.

## Running your first test

The fastest way to confirm the skill is loaded is to ask it a few short questions that do not require your own data. Try any of these:

- "Use the GAPIT AI skill to explain what Multi_iter does in GAPIT."
- "Use the GAPIT AI skill to generate a standard BLINK GWAS code template."
- "Use the GAPIT AI skill to explain the RefInf column in the BLINK prediction output."

A correct answer, in the form of an explanation or a code template, means the skill is active. If you get nothing, check that you opened the correct workspace folder, confirm the `SKILL.md` file exists, then restart TRAE and reopen the workspace.

To run a real analysis, ask the skill to build a GWAS or GS pipeline for your data. It will produce an R script based on one of its mode-code templates. Save the script, point it at your own genotype and phenotype files, and run it in R or RStudio. A minimal check that GAPIT loads correctly is to source the bundled function file and confirm the main function exists:

```r
skill_dir <- "path/to/your/unzipped/folder"
source(file.path(skill_dir, "data/gapit-resources/code/gapit_functions.txt"))
exists("GAPIT")   # should print TRUE
```

If that prints `TRUE`, your environment is ready for analysis.

## The two mode-code templates

The package ships two reference templates under `data/gapit-resources/code/GAPIT AI coding/`. `GAPIT_AI_GWAS_template.R` is the complete GWAS pipeline, and `GAPIT_AI_GS_template.R` is the complete GS pipeline. The skill treats these as the canonical structure for any code it generates.

The GWAS template walks through quality control, population structure, optional phenotype simulation, model fitting under GLM, MLM, CMLM, MLMM, SUPER, FarmCPU, and BLINK, extraction of significant markers, and power and false-discovery-rate calculation using GAPIT's own `GAPIT.FDR.TypeI` function.

The GS template walks through quality control, optional simulation, repeated k-fold cross-validation across gBLUP, cBLUP, sBLUP, GAGBLUP, MAS, BayesB, and rrBLUP, and reports two accuracy measures: correlation with the observed phenotype and correlation with the true breeding value.

Both templates run a lightweight quality control pass before calling GAPIT and print a summary of what was filtered, so you always know the analysis ran on quality-controlled data. The templates are a starting point, not a locked contract. You are free to edit them to match your files and study design.

## How the R code loads GAPIT

Generated scripts never call `library(GAPIT)`. They load GAPIT from the function file that ships with the package:

```r
source(file.path(skill_dir, "data/gapit-resources/code/gapit_functions.txt"))
```

This is intentional. The bundled file removes the dependency on a separately installed GAPIT R package and avoids version mismatches between your local copy and the skill. When you receive a script, keep that source line at the top and only change the file paths and parameters below it.

## Estimating token cost

TRAE charges by token usage. Every model request is metered on its own, and its cost is computed as the number of tokens multiplied by that model's list price, which is quoted per one million tokens. The resulting amount is deducted from your plan's monthly Dollar Usage allowance.

Input tokens, output tokens, and cached tokens are priced separately. Output tokens generally cost several times more than input tokens, and reading from the model's cache is far cheaper than fresh input. A long conversation that repeatedly resends a large file as context therefore consumes much more than a short, focused prompt.

TRAE offers five plans. The Free plan provides a limited monthly allowance and only Auto mode. The paid plans, Lite, Pro, Pro+, and Ultra, add monthly Dollar Usage allowances of 8, 20, 60, and 200 US dollars respectively, along with manual model selection and faster queues. When your allowance runs out, you can enable On-Demand Usage to keep going and pay for actual token consumption in small increments. Prices and model rates change over time, so confirm the current numbers on the pricing page at `https://www.trae.ai/pricing` and the model pricing list at `https://docs.trae.ai/ide/models`.

Three habits keep GAPIT AI affordable. Use a smaller or cheaper model for routine explanation and reserve the large models for genuinely hard reasoning. Ask for a script once and run it locally in R, rather than asking the assistant to rerun the same analysis repeatedly in chat. And keep your inputs lean: give the assistant a file path instead of pasting an entire genotype matrix into the conversation.

## API pathway: bringing your own model

If you already pay for a model API from another provider, you can point TRAE at it and spend that quota instead of your TRAE allowance. Open Settings, then Models, then Add Model. Two paths are available.

The first path uses a preset provider. Choose one of TRAE's built-in providers, pick a model, and paste your API key. The second path uses a custom configuration. Choose the API format, which is either OpenAI Chat Completions or Anthropic Messages, then enter the request URL, the model ID, and your API key. Most OpenAI-compatible services, including DeepSeek, OpenRouter, and many proxies, work with the OpenAI format.

After adding a model, it appears in the model list at the bottom-right of the chat input, where you can select it for your GAPIT AI conversation. This is the cleanest way to keep using your own DeepSeek or similar quota while working inside TRAE.

## The feedback system

GAPIT AI improves through user feedback. When something goes wrong, whether it is an installation error, a confusing answer, or a wrong parameter, report it rather than working around it silently. The maintainers collect feedback through the GAPIT community group, GitHub issues, and email.

A useful report contains four pieces of information. State what kind of problem it is, such as installation, error, feature request, or knowledge gap. Describe concretely what happened. Say what you expected the assistant to do instead. And list your environment, meaning the TRAE version, the operating system, and the skill version.

The maintainers turn these reports into updated rules and question-bank entries in the next release, so a problem reported once stops recurring for everyone else.

## Getting help

The official TRAE documentation is at `https://docs.trae.ai/`. The models and custom-model guide is at `https://docs.trae.ai/ide/models`, and the plans and billing guide is at `https://docs.trae.ai/ide/new-plans-and-billing`. For GAPIT itself, refer to the project's GitHub repository and the papers it cites. For questions specific to this skill, contact the maintainers through the feedback channels above.
