---
layout: default
title: Home
---

# README (rendered with MathJax)

> This page mirrors your README so formulas render on the site.

{% raw %}
{% capture readme %}{% include_relative ../README.md %}{% endcapture %}
{{ readme }}
{% endraw %}
